**********************************************************************************************************************;
*** 1.Sampling ***;
**********************************************************************************************************************;
** Macro Sample_Core: Take a random sample **;
%macro Sample_Core(in=,out=,response_var=,samplesize=,ratio=,seed=);
	%local default_seed;
	%let default_seed=;
	%if %length(&seed) ne 0 %then %let default_seed=&seed;

	data _tmp0_;
		min_num_nonresp=int((&min_num_resp/&Ratio)-&min_num_resp);
		call symputx('default_num_response',&min_num_resp);
		call symputx('default_num_nonresponse',min_num_nonresp);
	run;

	data _tmp_;
		n1=int(&samplesize*&ratio);
		n2=&samplesize-n1;
		if n1 < &default_num_response then do;
			n1=&default_num_response;

			/* output warning message here*/
			put "Warning: Too few responses. Number of Responses is now set to &default_num_response..";
		end;
		if n2 < &default_num_nonresponse then do;
			n2=&default_num_nonresponse;
			put "Warning: Too few non-responses. Number of Non-Responses is now set to &default_num_response..";
		end;
		call symputx('_num_response_sample_',n1);
		call symputx('_num_nonresponse_sample_',n2);
	run;

	%put &_num_response_sample_ &_num_nonresponse_sample_;

	data _in1_ (compress=yes  pointobs=yes reuse=no) _in2_ (compress=yes  pointobs=yes reuse=no);
		set &in;
		if &response_var>0 then output _in1_;
			else output _in2_;
	run;

	data _null_;
		set _in1_ nobs=n1;
		set _in2_ nobs=n2;
		call symputx('_num_response_original_',n1);
		call symputx('_num_nonresponse_original_',n2);
		stop;
	run;

	%put &_num_response_original_ &_num_nonresponse_original_;

	%* Copy the same data;
	%if (&_num_response_original_=&_num_response_sample_) %then %do;
		data _out1_;
			set _in1_;
		run;
	%end;
	%else %if (&_num_response_original_>&_num_response_sample_) %then /* Sample */
		%Sample_Option_1(in=_in1_,out=_out1_,samplesize=&_num_response_sample_,seed=&default_seed);
	%else /* Bootstrap */
		%Sample_Option_2(in=_in1_,out=_out1_,samplesize=&_num_response_sample_,seed=&default_seed);

	%if (&_num_nonresponse_original_=&_num_nonresponse_sample_) %then %do;
		data _out2_;
			set _in2_;
		run;
	%end;
	%else %if (&_num_nonresponse_original_>&_num_nonresponse_sample_) %then /* Sample */
		%Sample_Option_1(in=_in2_,out=_out2_,samplesize=&_num_nonresponse_sample_,seed=&default_seed);
	%else /* Bootstrap */
		%Sample_Option_2(in=_in2_,out=_out2_,samplesize=&_num_nonresponse_sample_,seed=&default_seed);

	data &out;
		set _out1_ _out2_;
	run;

%mend Sample_Core;

** If the size of desired sample is smaller than the size of the input dataset **;
%macro Sample_Option_1(in=,out=,samplesize=,seed=);
	%* if seed not specified;
	%if %length(&seed)=0 %then %do;
		proc surveyselect data=&in out=&out n=&samplesize noprint;
		run;
	%end;
	%else %do;
		proc surveyselect data=&in out=&out n=&samplesize seed=&seed noprint;
		run;
	%end;
%mend Sample_Option_1;

** If the size of the desired sample is greater then the size of input **;
%macro Sample_Option_2(in=,out=,samplesize=,seed=);
	data _NULL_;
		set &in nobs=_nsize_;
		Extra_num=&samplesize-_nsize_;
		call symputx('Extra_num',Extra_num);
	run;

	data &out;
		set &in nobs=_nsize_;
		do _i_=1 to &Extra_num;
			_tmp_index_=ceil(_nsize_*ranuni(&default_seed));
			set &in point=_tmp_index_;
			output;
		end;
		drop _i_;
		stop;
	run;

	data &out;
		set &in &out;
	run;

%mend Sample_Option_2;

** Master Sampling macro **;
%macro CE_Sampling(inds=, outds=OUT.CE1_Resampled);
	%*Split the original datasets;
	data _mod _test;
		set &inds;
		%if  &exclusion_if NE  %then %do;
			&exclusion_if then delete;
		%end;
		&split_if then mod_val_test=1;
			else mod_val_test=3;
		%if %upcase(&ds_present)=Y %then %do;
			%inc "&path_DS/DS Recodes.txt";
		%end;
		if mod_val_test=3 then output _test;
			else output _mod;
	run;

	%*bootstrap the modeling portion if requested;
	%if %lowcase(&Binary_dv)=y and %lowcase(&bootstrap)=y %then %do;

		%*Get key components for next sampling macro;
		proc means data=_mod(keep=&dep_var) n mean noprint;
			var &dep_var;
			output out=_temp  n=count1 mean=true_resp;
		run;

		data _NULL_;
			set _temp;
			num_resp=count1*true_resp;
			if num_resp>=&min_num_resp then do;
				spsize=int((count1*true_resp)/&oversampled_rr);
			end;
			else if num_resp<&min_num_resp then do;
				spsize=int(&min_num_resp/&oversampled_rr);
				put "Warning: Too few responses. Number of Responses is now set to &min_num_resp...";
			end;
			call symputx('spsize',spsize);
			call symputx('true_resp',true_resp);
		run;

		%put ** New modeling sample size: &spsize;
		%put ** Oversampled Response Rate: &oversampled_rr;
		%put ** Original Resp Rate: &true_resp;

		%Sample_Core(in=_mod,out=_mod,response_var=&dep_var,samplesize=&spsize,ratio=&oversampled_rr,seed=&seed);

		%*QC;
		proc freq data=&inds;
			title "&dep_var Frequency in Original Dataset";
			tables &dep_var;
		run;

		proc freq data=_mod;
			title "&dep_var Frequency in Oversampled Training Dataset";
			tables &dep_var;
		run;

		proc freq data=_test;
			title "&dep_var Frequency in unsampled Validation Dataset";
			tables &dep_var;
		run;
		title ' ';
	%end;

	%*combine splitted/oversampled data;
	data &outds;
		set _mod  _test;
		if mod_val_test=1 and ranuni(9545)>0.5 then
			mod_val_test=2;
	run;

	%*Output the sample response rate;
	proc means data=&outds nway;
		class mod_val_test;
		var &dep_var;
		output out=out.CE1_Sample_Rate mean=;
	run;

	data out.CE1_Sample_Rate;
		length partition $ 5;
		set out.CE1_Sample_Rate (drop=_type_);
		if mod_val_test=1 then partition='Model';
			else if mod_val_test=2 then partition='Val';
			else if mod_val_test=3 then partition='Test';
	run;

	%*clean up space;
	proc datasets lib=work nolist;
		delete _tmp_ _in1_ _in2_ _out1_ _out2_  _mod _test;
	quit;

%mend;