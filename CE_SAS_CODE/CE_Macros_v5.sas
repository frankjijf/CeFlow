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

**********************************************************************************************************************;
*** 2.EDA,profiling and Recode ***;
**********************************************************************************************************************;
** Common ordinal and continuous processing macro **;
%macro pnum(insdn,var,nmiss,typ);
    %* Get Non-missing key statistics;
    proc means data=&insdn mean median min p1 p25 p75 p99 max NOPRINT;
      var &var;
      output out=EDA mean=var_mean median=var_median min=var_min p1=var_p1 p25=var_p25 p75=var_p75 p99=var_p99 max=var_max / noinherit;
    run;
	%let skip = 1;
    data EDA; set EDA;
      %* calculate upper and lower bounds;
	  iqr=Max((var_p75-var_p25),(var_p99-var_p75),(var_p25-var_p1));
	  var_lb=Min(Max((var_p25- 1.5*iqr),var_min),var_p1);
	  var_ub=Max(Min((var_p75+ 1.5*iqr),var_max),var_p99);
      if var_lb=var_ub then do;
        var_lb=var_min;
        var_ub=var_max;
	  end;
	  var_mid = (var_max - var_min)/2;

      CALL SYMPUTX('var_LB' ,var_lb);
      CALL SYMPUTX('var_UB' ,var_ub);
      CALL SYMPUTX('var_median' ,var_median);

	  if upcase("&&impmethod&typ") in ('MEAN','STD') then CALL SYMPUTX('var_miss' ,var_mean);
	    else if upcase("&&impmethod&typ") in ('MEDIAN','IQR','MAD') then CALL SYMPUTX('var_miss' ,var_median);
	    else if upcase("&&impmethod&typ") in ('RANGE') then CALL SYMPUTX('var_miss' ,var_min);
	    else if upcase("&&impmethod&typ") in ('MIDRANGE') then CALL SYMPUTX('var_miss' ,var_mid);
	    else if upcase("&&impmethod&typ") in ('SUM','EUCLEN','USTD','MAXABS') then CALL SYMPUTX('var_miss' ,0);
	    else CALL SYMPUTX('skip',0);
	  
    run;

	%* Cap_FLR before transformation;
    data NONMISSING; set &insdn (where=(missing(&var)=0));
      &var=MIN(MAX(&var,&var_lb),&var_ub);
	  %if &&transformation&typ=Y %then %do;
        SQ_&var=&var**2;               
        SR_&var=sqrt(max(&var,0));     
        LN_&var=log(max(&var,0.00001)); 
	    *IV_&var= -1/max(&var,0.00001);    
        *EP_&var= -exp(Min(0,-&var));
        keep &dep_var &var SQ_&var SR_&var LN_&var /*IV_&var EP_&var*/;
		%let mod_list=&var SQ_&var SR_&var LN_&var;
      %end;
      %else %do;
        keep &dep_var &var;
		%let mod_list=&var;
	  %end;
	run;

	%let ck = %sysevalf(&nmiss - &min_size < 0);
	%if &skip = 0 and %upcase(&&impmethod&typ) ^= ER %then %do;
	  proc stdize data=&insdn out=dummy outstat=tmp method = &&impmethod&typ;
        var &var;
      run;
	  data _NULL_; set tmp (where=(_type_="LOCATION"));
	    call symputx("var_miss",&var);
	  run;
	  proc datasets nolist; delete tmp dummy; quit; run;
	%end;
	%else %if %upcase(&&impmethod&typ) = ER and &ck = 0 %then %do;
	  %* Get missing data average target;
      proc means data=&insdn mean noprint;
        var &dep_var;
		where missing(&var);
        output out=M_RR mean=Missing_RR;
      run;
	  data _NULL_; set M_RR;
	    %if %upcase(&Binary_dv) = Y %then %do;
          if Missing_RR=0 then Missing_RR=0.0001;
            else if Missing_RR=1 then Missing_RR=0.9999;
		%end;
        call symputx('Missing_RR' ,Missing_RR);
      run;
	  proc datasets nolist; delete M_RR; quit; run;
	%end;

	%if %upcase(&Binary_dv) = Y %then %do; 
      %* Run univariate logistic regression;
	  ods listing close;                                                                             
	  ods output parameterestimates=parm association=assoc;                                          
	  proc logistic data=NONMISSING desc namelen=32;
 		model &dep_var=&mod_list /PARMLABEL selection=forward STOP=1 slentry=1;
      run;
	  ods listing; 
	  %if %sysfunc(exist(assoc)) %then %do; 
	    data _NULL_; set assoc(keep=cvalue1 cvalue2);
	      if _n_ = 1 then call symputx("Concordant",cvalue1);
	        else if _N_ = 4 then call symputx("CValue",cvalue2);
	    run;
		proc datasets nolist; delete assoc; quit; run;
	  %end;
	  %else %do;
	    %let Concordant = . ;
		%let CValue = . ;
	  %end;
	%end;
	%else %do;
	  %* Run univariate simple regression;
	  ods listing close;
      ods output SelParmEst=parm SelectionSummary=SelectionSummary;
	  proc reg data=NONMISSING;
	    model &dep_var=&mod_list/selection=forward MAXSTEP=1 slentry=0.999;
	  run;
 	  ods listing;
	  %if %sysfunc(exist(SelectionSummary)) %then %do; 
	    data _NULL_; set SelectionSummary(keep=ModelRsquare obs=1);
	      call symputx("RSquare",ModelRsquare);
	    run;
		proc datasets nolist; delete SelectionSummary; quit; run;
	  %end;
	  %else %do;
	    %let RSquare = . ;
	  %end;
	%end;
 
    data _NULL_; set parm(obs=1);
	  call symputx('Intercept',Estimate);
    run; 

	 data _NULL_; set parm (firstobs=2 obs=2);
	 	%if %upcase(&Binary_dv) = Y %then %do;
	 		call symputx('prob',ProbChiSq);
		%end;
		%else %do;
	 		call symputx('prob',ProbF);
		%end;
	 run;
	 %let ckP = %sysevalf(&prob - .05 > 0);

	data parm; set parm (firstobs=2 obs=2);
	  length relationship $ 6  ck 3.; 

	  _Trans=substr(variable,1,3);
	  if _Trans not in ('SQ_','SR_','LN_','IV_','EP_') or (variable = "&var") then _Trans='';

      %if %upcase(&&impmethod&typ) = ER %then %do;
	    %if %upcase(&Binary_dv) = Y %then %do;
	      %*missing imputation based on missing target rate;
	      %if &ck = 1 or &ckP = 1 %then %do;
				miss_impute=&var_median;
			%end;
		    %else %do;
		           if _Trans='SQ_' then miss_impute=SQRT(Max((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate,0));
 		      else if _Trans='SR_' then miss_impute=((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate)**2;
		      else if _Trans='LN_' then miss_impute=exp((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate);
		      else if _Trans='IV_' then miss_impute=-1/((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate);
		      else if _Trans='EP_' then miss_impute=-log(Max(-(log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate,0.00001));
              else                      miss_impute=(log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate;
          %end;
        %end;
	    %else %do;
	      %*missing imputation based on missing target rate;
		   %if &ck = 1 or &ckP = 1 %then %do;
				miss_impute=&var_median;
			%end;
  		   %else %do;
		           if _Trans='SQ_' then miss_impute= SQRT(Max((&Missing_RR-&Intercept)/Estimate,0));
		      else if _Trans='SR_' then miss_impute= ((&Missing_RR-&Intercept)/Estimate)**2;
		      else if _Trans='LN_' then miss_impute= exp((&Missing_RR-&Intercept)/Estimate);
		      else if _Trans='IV_' then miss_impute= -1/((&Missing_RR-&Intercept)/Estimate);
			   else if _Trans='EP_' then miss_impute= -log(Max(-(&Missing_RR-&Intercept)/Estimate,0.00001));
              else                      miss_impute= (&Missing_RR-&Intercept)/Estimate;
         %end;
		%end;
		if miss_impute<&var_LB then miss_impute=&var_LB;
		  else if miss_impute>&var_UB then  miss_impute=&var_UB;
		call symputx("var_miss",miss_impute);
      %end;

      if Estimate>0 then sign='(+)'; else sign='(-)';
      relationship=compress(_Trans||sign,'');

	  ** Output fields to append to the vars table;
	  call symputx("New_var",variable);
	  call symputx("Sign",sign);
	  call symputx("Relation",relationship);
	  %if %upcase(&Binary_dv) = Y %then %do;
	    call symputx("Prob",ProbChiSq);
	    call symputx("Chisq",WaldChiSq);
	  %end;
	  %else %do;
	    call symputx("Prob",ProbF);
		call symputx("FValue",FValue);
	  %end;
	run;           
     
	%** Update vars table **;
	data vars2; set vars2;
	  if name = "&var" then do;
        var_lb = &var_lb;
        var_ub = &var_ub;
        var_median = &var_median;
		miss_impute = &var_miss;
		%if %upcase(&Binary_dv) = Y %then %do;
		  Concordant = &Concordant;
		  CValue = &CValue;
		  Chisq = &Chisq;
		  PValue = &Prob;
		%end;
		%else %do;
		  RSquare = &RSquare;
		  FValue = &FValue;
		  PValue = &Prob;
		%end;
		if PValue <= &PValue then do;
		  new_var = "&Prefix"||"&New_var";
		  Sign = "&Sign";
		  Relationship = "&Relation";
		end;
      end; 
	run;

	%if %sysevalf(&Pvalue-&prob>=0) %then %do;
	  %if %upcase("&&stdmethod&typ") ^= "NO" %then %do;
	    %* Standardization;
	    data temp; set &insdn (keep=&var);
	      _Trans=substr("&Relation",1,3);
	      &New_var = &var;
		  if missing(&New_var) then &New_var = &var_miss;
		  %if %upcase(&&cap_flr&typ) = Y %then %do;
	        &New_var = MIN(MAX(&New_var, &Var_LB), &var_UB);
	      %end;
		  *standard transformation;
	           if _Trans='SQ_' then &New_var = &New_var **2;
          else if _Trans='SR_' then &New_var = SQRT(MAX(&New_var,0));
	      else if _Trans='LN_' then &New_var = LOG(MAX(&New_var,0.00001));
	      else if _Trans='IV_' then &New_var = -1/(MAX(&New_var,0.00001));
	      else if _Trans='EP_' then &New_var = -EXP(MIN(-&New_var,0));
	    run;
		proc stdize data=temp out=dummy outstat=std_tmp method = &&stdmethod&typ;
          var &New_var;
        run;
		data _NULL_; set std_tmp;
		  if _type_ = 'LOCATION' then call symputx("loc",&new_var);
		  if _type_ = 'SCALE' then call symputx("scale",&new_var);
		run;
		%** Update vars table **;
	    data vars2; set vars2;
	      if name = "&var" then do;
            Location = &loc;
            Scale = &scale;
		  end;
	    run;
		proc datasets nolist; delete temp dummy std_tmp; quit; run;
	  %end;
	  data _NULL_; set vars2 (where=(name="&var"));
	    _Trans=substr(relationship,1,3);
		if missing(label) then lab2 = name;
		  else lab2 = label;

		%if &typ = O %then %do;
          FILE  "&Path_output.CE2_Ordinal_Var_Recode.txt" mod;
		%end;
		%else %do;
          FILE  "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
		%end;
	    PUT ' ' '0d'x;
	    PUT "*** RECODE " name ": " label " ***;" '0d'x;
	    * Missing imputation;
        Put "IF missing(" name +(-1) ") THEN  &PREFIX" name " = " Miss_Impute " ; " '0d'x;
		* Valid values;
		PUT @2 "ELSE &PREFIX" name " = " name ";" '0d'x;
	    * Capping and flooring;
	    %if %upcase(&&cap_flr&typ) = Y %then %do;
	      Put "&PREFIX" name " = MIN(MAX(&PREFIX" name ", " Var_LB +(-1) "), " var_UB +(-1) ");" '0d'x;
	    %end;
        * Untransformed variable label;
		%if &typ = O %then %do;
	      PUT @2  "Label &PREFIX" name ' = "' lab2 ': Ordinal Recode ' Sign '";' '0d'x;
		%end;
		%else %do;
		  PUT @2  "Label &PREFIX" name ' = "' lab2 ': Continuous Recode ' Sign '";' '0d'x;
		%end;

	    *standard transformation;
	         if _Trans='SQ_' then do; Put new_var "= &PREFIX" name "**2;" '0d'x;  
                                      Put @2 " Label " new_var '= "' lab2 'SQUARE ' Sign '";' '0d'x; end;
        else if _Trans='SR_' then do; Put new_var " = SQRT(MAX(&PREFIX" name ",0));" '0d'x;   
												  Put @2 " Label " new_var '= "' lab2 'SQRT ' Sign '";' '0d'x; end;
	    else if _Trans='LN_' then do; Put new_var "= LOG(MAX(&PREFIX" name ",0.00001));" '0d'x;   
												  Put @2 " Label " new_var '= "' lab2 'LOG ' Sign '";' '0d'x; end;
	    else if _Trans='IV_' then do; Put new_var "= -1/(MAX(&PREFIX" name ",0.00001));" '0d'x;   
												  Put @2 " Label " new_var '= "' lab2 'NEGATIVE INVERSE ' Sign '";' '0d'x; end;
	    else if _Trans='EP_' then do; Put new_var "= -EXP(MIN(-&PREFIX" name ",0));" '0d'x;   
												  Put @2 " Label " new_var '= "' lab2 'EXPONENTIAL ' Sign '";' '0d'x; end;
		* Standardization;
	    %if %upcase("&&stdmethod&typ") ^= "NO" %then %do;
		  if (location ^= . and scale ^= .) then put new_var "= (" new_var "- (" location +(-1) ") ) / " scale ";" '0d'x;       
	    %end;

	  run;
 	%end;
	proc datasets nolist; delete parm nonmissing; quit; run;
%mend pnum;
** Ordinal and continuous profiling macro 1 **;
%macro prof1(insdn,var);
  proc summary data=&insdn nway missing;
    var &dep_var;
	class &var;
	output out=prof (drop=_type_ rename=_freq_=xcount) mean=xmean;
  run;
  data prof; set prof;
    length xcategory $256.;
	if missing(&var) then xcategory = 'Missing';
	  else xcategory = trim(left(put(&var,best8.)));
  run;
%mend;
** Ordinal and continuous profiling macro 2 **;
%macro prof2(insdn,var);
  %if %upcase(&equal_dist) = Y %then %do;
    data _NULL_; set EDA;
	  range = (var_p99 - var_p1)/&num_category;
	  call symputx('cut_lo',var_p1);
	  call symputx('cut_hi',var_p99);
	  call symputx('range',range);
    run;
	%put &range;
	data tmp; set &insdn (keep=&dep_var &var);
	  if missing(&var) then bin = .;
	    else if &var < &cut_lo then bin = 1;
		else if &var >= &cut_hi then bin = input("&num_category",best10.);
		else do;
		  do k = 1 to &num_category;
		    if &var>= &cut_lo+(k-1)*&range and &var< &cut_lo+k*&range then bin=k;
		  end;
		end;
	run;
  %end;
  %else %do;
    proc rank data=&insdn (keep=&dep_var &var) out=tmp ties=High group=&num_category;
      var &var;
	  ranks bin;
    run;
  %end;
  proc summary data=tmp nway missing;
    var &var &dep_var;
	class bin;
	output out=prof (drop=_type_ rename=_freq_=xcount) min(&var)=lo max(&var)=hi mean(&dep_var)=xmean;
  run;
  data prof (drop=tag); set prof end=eof;
    length xcategory $256.;
	retain tag 0;
	if missing(bin) then do;
      xcategory = 'Missing';
	  tag = 1;
	end;
	else if _N_ = 1 or (_N_ = 2 and tag = 1) then xcategory = "Low to " || trim(left(put(hi,best8.)));
	else if eof then xcategory = trim(left(put(lo,best8.))) || " to High";
	else xcategory = trim(left(put(lo,best8.))) || " to " || trim(left(put(hi,best8.)));
  run;
  proc datasets nolist; delete tmp eda; run;
%mend;
** Ordinal and continuous profiling macro 3 **;
%macro prof3(typ);
  proc sql;
    create table prof2 as
	select a.name as variable, a.label, b.xcategory,
	       b.xcount, b.xmean
    from vars2 a, prof b
	where a.name = "&var";
  quit;
  data prof3 (drop=xcount xmean xcategory); set prof2 end=eof;
    length category $256.;
    length type $10.;
	length star $8.;
	if "&typ" = "C" then type = 'Continuous';
	  else type = 'Ordinal';
    if _N_ = 1 then do;
      category = "Overall";
      Average_DV = &overall_avg;
      Count = &nobs;
      Percent = 1;
      index = 100;
      output;
    end;
	category = xcategory;
	count = xcount;
	percent = xcount / &nobs;
	Average_DV = xmean;
    index = (Average_DV / &overall_avg)*100;
    if index >= 110 then star = '* (+)';
      else if index > 100 then star = '  (+)';
      else if index <= 90 then star = '* (-)';
      else if index <= 100 then star = '  (-)';
      else star = '  (0)';
	  output;
	run;
	data out.CE2_profile; set out.CE2_profile prof3; run;
	proc datasets nolist; delete prof prof2 prof3; run;
%mend;

** Binary variable processing macro **;
%macro pbin(insdn,var,typ);
  %** Get missing count **;
  proc sql noprint;
    select max(sum(missing(&var)),0) into : miss_ck from &insdn;
  quit;
  %put processing binary variable &var with type &typ;
  %if %sysevalf((&nobs*&missrate)>=&miss_ck) %then %put missing count is &miss_ck;
  %else %put missing count is &miss_ck which is too high;
  %** If missing count is acceptable **;
  %if %sysevalf((&nobs*&missrate)>=&miss_ck) %then %do;
    %** Get counts of values that should be coded as 1 **;
    %if &typ = 1 %then %do;  ** Numeric variables **;
      proc sql noprint;
        select sum(case when &var = 1 then 1 else 0 end) into : cnt_ck from &insdn;
      quit;
	%end;
	%else %do;  %** Character variables **;
      proc sql noprint;
        select sum(case when strip(&var) in ('1','Y','y') then 1 else 0 end) into : cnt_ck from &insdn;
      quit;
	%end;
    %if %sysevalf((&nobs*(1-&concrate)) <= &cnt_ck and &cnt_ck <= (&nobs*&concrate)) %then %do;
	  %** If count is okay, write out code to file **;
	  data _null_; set vars (where=(name="&var"));
		FILE  "&Path_output.CE2_Binary_Var_Recode.txt" mod;
		PUT ' ' '0d'x;
		PUT "*** RECODE " name ": " label " ***;" '0d'x;
		%if &typ = 1 %then %do;
          PUT "if " name  " = 1 then &PREFIX" name " = 1; else &PREFIX" name " = 0;" '0d'x;
		%end;
		%else %do;
          PUT "if " name  " in ('1','Y','y') then &PREFIX" name " = 1; else &PREFIX" name " = 0;" '0d'x;
		%end;
		if missing(label) then PUT @2  "Label &PREFIX" name ' = "' name ': Binary Recode";' '0d'x;
		  else PUT @2  "Label &PREFIX" name ' = "' label ': Binary Recode";' '0d'x;
	  run;
      %if %upcase(&profiling) = Y %then %do;
	    %if &typ = 1 %then %do;
  	      data tmp1; set &insdn (keep=&dep_var &var);
		    if &var = 1 then newvar = 1;
		      else newvar = 0;
		  run;
		%end;
		%else %do;
		  data tmp1; set &insdn;
		    if &var in ('1','Y','y') then newvar = 1;
		      else newvar = 0;
		  run;
		%end;
		proc sql;
		  create table tmp2 as
		  select a.name as variable, a.label, b.newvar,
		         count(*) as xcount,
		         mean(b.&dep_var) as xmean
		  from vars a, tmp1 b
		  where a.name = "&var"
		  group by a.name, a.label, b.newvar
		  order by a.name, a.label, b.newvar;
		quit;
	    data tmp3 (drop=xcount xmean newvar); set tmp2 end=eof;
		  by newvar;
		  length type $10.;
	      length category $256.;
	  	  length star $8.;
		  type = 'Binary';
          if _N_ = 1 then do;
            category = "Overall";
            Average_DV = &overall_avg;
            Count = &nobs;
            Percent = 1;
            index = 100;
            output;
          end;
		  count = xcount;
		  percent = xcount / &nobs;
		  Average_DV = xmean;
		  index = (Average_DV / &overall_avg)*100;
		  %if &typ = 1 %then %do;
            if newvar=0 then category = "Missing,0";
          	  else category = "1";
		  %end;
		  %else %do;
		    if newvar=0 then category = "Missing,0,N";
          	  else category = "1,Y";
		  %end;
          if index >= 110 then star = '* (+)';
            else if index > 100 then star = '  (+)';
            else if index <= 90 then star = '* (-)';
          	else if index <= 100 then star = '  (-)';
          	else star = '  (0)';
		    output;
	     run;
		 data out.CE2_profile; set out.CE2_profile tmp3; run;
		 proc datasets nolist; delete tmp: ; run;
      %end;
	%end;
	%else %do;  %** If variable is too concentrated, write note to log **;
	  %put variable &var is too concentrated to use;
	%end;
	%** Add counts to summary file **;
	data vars; set vars;
      if name = "&var" then do;
	    miss_cnt = &miss_ck;
	    ones_cnt = &cnt_ck;
		%if %sysevalf((&nobs*(1-&concrate)) <= &cnt_ck and &cnt_ck <= (&nobs*&concrate)) %then %do;
          good = 1;
		%end;
	  end;
    run;
  %end;
  %** If missing count is too high **;
  %else %do;
    data vars; set vars;
      if name = "&var" then do;
	    miss_cnt = &miss_ck;
  	  end;
    run;
  %end;
  
%mend pbin;
 
** Binary processing control macro **;
%macro bin_cntl(insdn);
%** Initialize code file **;
  data _null_;
    FILE  "&Path_output.CE2_Binary_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                        CE2_BINARY_VAR_RECODE.TXT                        ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
  run;
%** Get variables to process **;
  proc contents data=&insdn (keep=&binvar)
      out=vars (keep=varnum name type length label)
      noprint;
  run;
  proc sql noprint;
    select count(*) into : ck from vars;
  quit;
  %put number of binary variables = &ck;

%** If there are binary variables to process **;
  %if &ck > 0 %then %do;
    %** Add new variables to summary file **;
    data vars; set vars;
	  length miss_cnt ones_cnt good 8.;
	run;
	%** Create macro variables to control processing **;
    data _null_; set vars (where=(length(strip(name)) + length(strip("&prefix")) <= 32)) end=eof;
      if eof then call symputx("VARCNT",_N_);
      call symputx("IV"|| trim(left(put(_N_,4.)))  ,name);
	  call symputx("typ"|| trim(left(put(_N_,4.)))  ,type);
    run;
	%** Loop through for each variable **;
    %do _I_ = 1 %to &varcnt;
      %pbin(&insdn,&&iv&_I_,&&typ&_I_);
    %end;
	%** Create new variable name for surviving variables **;
	data out.CE2_binary_vars (drop=good); set vars;
	  length new_var $32.;
	  if good = 1 then new_var = "&prefix"||name;
	run;
	%** Write out list of variables to keep to code file **;
	%let rck = 0;
	data _NULL_; set out.CE2_binary_vars (where=(missing(new_var)=0)) end=eof;
	  m = mod(_N_,7);
      length v $256.;
      retain v ;
      if m = 1 then v = new_var;
        else v = strip(v) || " " || strip(new_var);
	  FILE  "&Path_output.CE2_Binary_Var_Recode.txt" mod;
	  if _N_=1 then do;
        PUT ' ' '0d'x;
		  PUT ' ' '0d'x;
		PUT '%LET KEEPLIST_B =';
      end;
      if m = 0 or eof then PUT v '0d'x;
	  if eof then PUT ';';
	  if eof then call symputx("rck" ,_N_);
    run;
	%put number of recoded binary variables = &rck;

	%* Capture old and new variable name;
	data tmp (keep=orig_var orig_label variable);
	  length variable orig_var $32.;
	  length orig_label $256.;
      set out.CE2_binary_vars (where=(missing(variable)=0) rename=(name=orig_var label=orig_label new_var=variable));
	run;
	data keep_vars; set keep_vars tmp; run;

	proc datasets nolist; delete vars tmp; quit; run;
  %end;
%** If there are no binary variables to process **;
  %else %do;
    data _null_;
      FILE  "&Path_output.CE2_Binary_Var_Recode.txt" mod;
      PUT ' ' '0d'x;
      PUT '%LET KEEPLIST_B = ;' '0d'x;
    run;
  %end;
%mend bin_cntl;

** Nominal variable processing macro **;
%macro pnom(insdn,var,typ);
  %** Summarize **;
  proc sql;
    create table tmp as
	select &var, count(*) as dcount, mean(&dep_var) as dmean, var(&dep_var) as dvar
	from &insdn
	group by &var
	order by dmean;
  quit;
  proc sql noprint;
    select sum(case when missing(&var) then dcount else 0 end) into : misscnt from tmp;
	select max(dcount) into : maxcnt from tmp;
	select count(*) into : unqcnt from tmp;
  quit;
  %** Test **;
  %put processing nominal variable &var with type &typ;
  %if %sysevalf((&nobs*&missrate)<&misscnt) %then %put missing count is &misscnt which is too high;
  %else %do;
    %put missing count is &misscnt;
    %if %sysevalf((&nobs*&concrate)<&maxcnt) %then %put variable &var is too concentrated to use;
    %else %do;
	  %if %sysevalf(&valcnt^=0 and &valcnt<=&unqcnt) %then %put variable &var has too many values to use;
	  %else %do;  %** Variable is usable;
        %** Collapse values based on counts;
        DATA tmp1; SET tmp END=eof;
          LENGTH tcount tmean tvar tgroup cumcnt tgrpcnt 8.;
          RETAIN tcount tmean tvar tgroup cumcnt tgrpcnt;
	      cumcnt + dcount;
          IF _n_ = 1 THEN DO;
            tcount = dcount;
			tmean = dmean;
			tvar = dvar;
            tgroup = 1;
		    tgrpcnt = 1;
          END;
          ELSE DO;
            IF tcount <= &minbinn or (cumcnt - dcount) >= (&nobs - &minbinn) THEN DO;
              tmean = ((tcount * tmean) + (dcount * dmean))/(tcount + dcount);
              tvar = (((tcount - 1) * tvar) + ((dcount - 1) * dvar))/(tcount + dcount - 2);
              tcount = tcount + dcount;
              tgroup = tgroup;
		      tgrpcnt = tgrpcnt + 1;
            END;
			ELSE DO;
			  tcount = dcount;
              tmean = dmean;
              tvar = dvar;
              tgroup = tgroup+1;
		      tgrpcnt = 1;
			END;
          END;
          KEEP &var dcount dmean dvar tgroup tcount tmean tvar tgrpcnt;
        RUN;
		%** Create dataset with final record for group only;
		data tmp2; set tmp1 (drop=&var dcount dmean dvar);
		  by tgroup;
		  if last.tgroup;
		run;

        %** Calculate bonferroni adjustment on alpha value;  
        %IF %UPCASE(&bonfer) = Y %THEN %DO;
          PROC SQL NOPRINT;
            SELECT COUNT(DISTINCT group)
            INTO :ncomps
            FROM tmp2;
          QUIT;
          %LET ncomps = %SYSEVALF(&ncomps - 1);
          %IF &ncomps > 1 %THEN %LET ftalpha = %SYSEVALF(&talpha./&ncomps);
            %ELSE %LET ftalpha  = &talpha;
        %END;
        %ELSE %LET ftalpha  = &talpha;

		%** Collapse values based on variance;
		DATA tmp3; SET tmp2 END=eof;
          LENGTH row fcount fmean fvar fgroup grpcnt 8.;
          RETAIN fcount fmean fvar fgroup grpcnt;
	      row = _N_;
          IF _n_ = 1 THEN DO;
            fcount = tcount;
            fmean = tmean;
            fvar = tvar;
            fgroup = 1;
		    grpcnt = tgrpcnt;
          END;
          ELSE DO;
            pvar = (((fcount - 1) * fvar) + ((tcount - 1) * tvar))/(fcount + tcount - 2);
            t = (fmean - tmean)/SQRT(pvar * ((1/fcount) + (1/tcount)));
            df = (fcount + tcount - 2);
            prob = (1-PROBT(abs(t),df));
            IF prob <= &ftalpha THEN DO;
              fcount = tcount;
              fmean = tmean;
              fvar = tvar;
              fgroup = fgroup+1;
		      grpcnt = tgrpcnt;
            END;
            ELSE DO;
              fmean = ((fcount * fmean) + (tcount * tmean))/(fcount + tcount);
              fvar = (((fcount - 1) * fvar) + ((tcount - 1) * tvar))/(fcount + tcount - 2);
              fcount = fcount + tcount;
              fgroup = fgroup;
		      grpcnt = grpcnt + tgrpcnt;
            END;
          END;
          KEEP row tgroup fgroup fcount fmean grpcnt;
        RUN;
        PROC SORT DATA = tmp3 out=tmp4; BY fgroup DESCENDING row; RUN;

        DATA tmp4 (drop=fcount fmean row grpcnt); SET tmp4;
          LENGTH xmean xcount group_size 8. name $32.;
          RETAIN xmean xcount group_size;
          BY fgroup DESCENDING row;
          IF FIRST.fgroup THEN DO;
            xmean = fmean;
            xcount = fcount;
		    group_size = grpcnt;
          END;
	      Index  = (xmean/&overall_Avg)*100 ;
	      diff = abs(xmean-&overall_Avg);
	      name="&var";
        RUN;

		%* Sort to keep smaller groups first;
	    PROC SORT DATA = tmp4; BY group_size DESCENDING diff xcount ; RUN;
        DATA _NULL_; SET tmp4 END=eof;
          IF eof THEN call symputx('last_group',fgroup);
        RUN;

		%* Combine data;
	    proc sql;
	      create table tmp5 as
	      select b.name, c.label, a.&var, b.fgroup, a.dmean, a.dcount, b.xmean,
                 b.xcount, b.diff, b.Index, b.group_size
	      from tmp1 a, tmp4 b, vars c
	      where a.tgroup = b.tgroup
	        and b.name = c.name
	      order by b.group_size, b.diff desc, b.xcount;
	    quit;

		%* Write out recode code;
		%if %upcase(&nom_method) = BINARY %then %do;
		  DATA _null_; SET tmp5 (where=(fgroup ^= &last_group)) END=eof;
            BY group_size DESCENDING diff ;
            FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" mod;

	        IF _n_ = 1 THEN DO;
		      PUT ' ' '0d'x;
		      PUT "*** RECODE " name ": " label "***;" '0d'x;
            END;
			%if &typ = 1 %then %do;
              IF FIRST.diff THEN PUT "&PREFIX" name +(-1) "_X" fgroup "= " name "in (" &var +(-1)@;  
              IF FIRST.diff NE 1 THEN PUT "," &var +(-1)@;
              IF LAST.diff THEN PUT ");" '0d'x;
			%end;
			%else %do;
              IF FIRST.diff THEN PUT "&PREFIX" name +(-1) "_X" fgroup "= " name "in ('" &var +(-1)@;  
              IF FIRST.diff NE 1 THEN PUT "','" &var +(-1)@;
              IF LAST.diff THEN PUT "');" '0d'x;
			%end;
          RUN;

          DATA _null_; SET tmp5 (where=(fgroup ^= &last_group)) END=eof;
            BY group_size DESCENDING diff ;
            FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" mod;

	        IF FIRST.diff THEN DO;
	          IF missing(label) then PUT @2  "Label &PREFIX" name +(-1) "_X" fgroup '= "' name ': Values ' &var +(-1)@;
		        ELSE PUT @2  "Label &PREFIX" name +(-1) "_X" fgroup '= "' label ': Values ' &var +(-1)@;
            END;
            IF FIRST.diff NE 1 THEN PUT "," &var +(-1)@;
            IF LAST.diff THEN PUT '";' '0d'x;
          RUN;
		  data nv (keep=name new_var); set tmp5 (where=(fgroup ^= &last_group));
		    by group_size DESCENDING diff ;
		    if first.diff;
			length new_var $32.;
			new_var = "&prefix" || strip(name) || "_X" || strip(left(put(fgroup,best2.)));
		  run;
		  data new_vars; set new_vars nv; run;
		%end;
		%else %do;
		  %if %upcase(&nom_method) = INDEX %then %let val =index;
		    %else %let val = xmean;
		  proc sort data=tmp5; by fgroup group_size &var; run;
		  proc sql noprint;
            select max(fgroup) into : last_group from tmp5;
          quit;
		  DATA _null_; SET tmp5 END=eof;
            BY fgroup group_size;
            FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" mod;

	        IF _n_ = 1 THEN DO;
	          IF fgroup = &last_group THEN STOP; 
		      ELSE DO;
		        PUT ' ' '0d'x;
		        PUT "*** RECODE " name ": " label "***;" '0d'x;
		      END;
            END;

	        IF fgroup = &last_group THEN DO;
			  IF first.group_size THEN PUT "else &PREFIX" name "= " &val ";" '0d'x;
			END;
	        ELSE DO;
			  %if &typ = 1 %then %do;
                IF _n_ = 1 AND FIRST.group_size THEN PUT "if " name "in (" &var +(-1)@;  
		          ELSE IF FIRST.group_size THEN PUT "else if " name "in (" &var +(-1)@;  
                IF FIRST.group_size NE 1 THEN PUT "," &var +(-1)@;
                IF LAST.group_size THEN PUT ") then &PREFIX" name "= " &val ";" '0d'x;
			  %end;
			  %else %do;
                IF _n_ = 1 AND FIRST.group_size THEN PUT "if " name "in ('" &var +(-1)@;  
		          ELSE IF FIRST.group_size THEN PUT "else if " name "in ('" &var +(-1)@;  
                IF FIRST.group_size NE 1 THEN PUT "','" &var +(-1)@;
                IF LAST.group_size THEN PUT "') then &PREFIX" name "= " &val ";" '0d'x;
			  %end;
	        END;
          RUN;

          DATA _null_; SET tmp5 (where=(fgroup ^= &last_group) obs=1) END=eof;
            FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" mod;
		    IF missing(label) then PUT @2  "Label &PREFIX" name ' = "' name ': Nominal Recode";' '0d'x;
		      ELSE PUT @2  "Label &PREFIX" name ' = "' label ': Nominal Recode";' '0d'x;
          RUN;
		  data nv (keep=name new_var); set tmp5 (where=(fgroup ^= &last_group) obs=1);
			length new_var $32.;
			new_var = "&prefix" || strip(name);
		  run;
		  data new_vars; set new_vars nv; run;
		%end;
		%if %upcase(&Profiling) = Y %then %do;
		  proc sort data=tmp5; by fgroup group_size &var; run;
		  data tmp6 (keep=type Variable label category count percent Average_DV index star); 
            set tmp5 (rename=(name=Variable)) end=eof;
            by fgroup;
			length type $10.;
            length category $256.;
            length star $8.;
            retain category;
			type = 'Nominal';
            if _N_ = 1 then do;
              category = "Overall";
          	  Average_DV = &overall_avg;
          	  Count = &nobs;
          	  Percent = 1;
          	  index = 100;
          	  output;
            end;
			count = xcount;
		    percent = xcount / &nobs;
		    Average_DV = xmean;
		    Index = (Average_DV / &overall_avg)*100;
            if first.fgroup and missing(&var) then category = "Missing";
              else if first.fgroup then category = strip(&var);
          	  else category = strip(category) || "," || strip(&var);
            if index >= 110 then star = '* (+)';
              else if index > 100 then star = '  (+)';
          	  else if index <= 90 then star = '* (-)';
          	  else if index <= 100 then star = '  (-)';
          	  else star = '  (0)';
            if last.fgroup then output;
          run;
		  data out.CE2_profile; set out.CE2_profile tmp6; run;
	    %end;
	  %end;
	%end;
  %end;
  %** Add counts to summary file **;
  data vars; set vars;
    if name = "&var" then do;
      unique_cnt = &unqcnt;
      miss_cnt = &misscnt;
      max_cnt = &maxcnt;
    end;
  run;
  proc datasets nolist; delete tmp: nv; quit; run;
 
%mend pnom;
 
** Nominal processing control macro **;
%macro nom_cntl(insdn);
%** Initialize code file **;
  data _null_;
    FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                        CE2_NOMINAL_VAR_RECODE.TXT                       ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
  run;
%** Get variables to process **;
  proc contents data=&insdn (keep=&nomvar)
      out=vars (keep=varnum name type length label)
      noprint;
  run;
  proc sql noprint;
    select count(*) into : ck from vars;
  quit;
  %put number of nominal variables = &ck;

%** If there are nominal variables to process **;
  %if &ck > 0 %then %do;
    %** Add new variables to summary file **;
    data vars; set vars;
	  length unique_cnt miss_cnt max_cnt 8.;
	run;

	%if %upcase(&nom_method) = BINARY %then %let maxlen=29;
	  %else %let maxlen=32;
    
	%** Create macro variables to control processing **;
    data _null_; set vars (where=(length(strip(name)) + length(strip("&prefix")) <= &maxlen)) end=eof;
      if eof then call symputx("VARCNT",_N_);
      call symputx("IV"|| trim(left(put(_N_,4.)))  ,name);
	  call symputx("typ"|| trim(left(put(_N_,4.)))  ,type);
    run;
	data new_vars; set _NULL_;
	  length name new_var $32.;
	run;
	%** Loop through for each variable **;
    %do _I_ = 1 %to &varcnt;
      %pnom(&insdn,&&iv&_I_,&&typ&_I_);
    %end;
	%** Create new variable name for surviving variables **;
	proc sort data=new_vars; by name; run;
    data nv; set new_vars;
      by name;
      length new_vars $256.;
      retain new_vars;
      if first.name then new_vars = new_var;
        else new_vars = strip(new_vars)||","||strip(new_var);
      if last.name then output;
    run;
	proc sql;
	  create table out.CE2_nominal_vars as
	  select a.*, b.new_vars
	  from vars a
	  left join nv b on a.name = b.name;
	quit;
	%** Write out list of variables to keep to code file **;
	%let rck = 0;
	data _NULL_; set new_vars (where=(missing(new_var)=0)) end=eof;
	  m = mod(_N_,7);
      length v $256.;
      retain v ;
      if m = 1 then v = new_var;
        else v = strip(v) || " " || strip(new_var);
	  FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" mod;
	  if _N_=1 then do;
        PUT ' ' '0d'x;
		  PUT ' ' '0d'x;
		PUT '%LET KEEPLIST_N =';
      end;
      if m = 0 or eof then PUT v '0d'x;
	  if eof then PUT ';';
	  if eof then call symputx("rck" ,_N_);
    run;
	%put number of recoded nominal variables = &rck;

	%* Capture old and new variable name;
	proc sql;
	  create table tmp as
	  select a.name as orig_var, a.label as orig_label, b.new_var as variable
	  from vars a, new_vars b
      where a.name = b.name;
	quit;
	data keep_vars; set keep_vars tmp; run;

	proc datasets nolist; delete vars new_vars nv tmp; quit; run;
  %end;
%** If there are no nominal variables to process **;
  %else %do;
    data _null_;
      FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" mod;
      PUT ' ' '0d'x;
      PUT '%LET KEEPLIST_N = ;' '0d'x;
    run;
  %end;
%mend nom_cntl;

** Ordinal variable processing macro **;
%macro pord(insdn,var,nmiss);
  %** Summarize **;
  proc sql;
    create table tmp as
	select &var, count(*) as count
	from &insdn
	group by &var;
  quit;
  proc sql noprint;
	select max(count) into : maxcnt from tmp;
	select count(*) into : unqcnt from tmp;
  quit;
  %* check if dep_var has the constant value in nonmissing part;
  proc sql noprint;
    select min(&dep_var)=max(&dep_var) into : cons_dep from &insdn
    where missing(&var)=0;
  quit;
  proc datasets nolist; delete tmp; quit; run;

  %* add to vars file;
  data vars2; set vars2;
    if name = "&var" then do;
	  unique_cnt = &unqcnt;
      max_cnt = &maxcnt;
    end; 
  run;

  %** Test **;
  %if %sysevalf((&nobs*&concrate)<&maxcnt) %then %put variable &var is too concentrated to use;
  %else %if &cons_dep = 1 %then %put &var has constant dependent value in nonmissing part and will be excluded;
  %else %do;
    %pnum(&insdn,&var,&nmiss,O);
	%** Profiling **;
    %if %upcase(&profiling) = Y %then %do;
      %if &unqcnt <= &num_category %then %do;
	    %prof1(&insdn,&var);
	  %end;
	  %else %do;
	    %prof2(&insdn,&var);
	  %end;
	  %prof3(O);
    %end;
  %end;

%mend pord;
 
** Ordinal processing control macro **;
%macro ord_cntl(insdn);
%** Initialize code file **;
  data _null_;
    FILE  "&Path_output.CE2_Ordinal_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                       CE2_ORDINAL_VAR_RECODE.TXT                        ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
  run;
%** Get variables to process **;
  proc contents data=&insdn (keep=&ordvar)
      out=vars (keep=varnum name type length label)
      noprint;
  run;
  proc sql noprint;
		select max(type) into : ck_type from vars;
  quit;
  %if &ck_type = 2 %then %do;
    	proc sql noprint; select name into : bad separated by ',' from vars where type = 2; quit;
		%put The following character variables were dropped from the ordinal variable list:;
		%put &bad;
		proc sql noprint; 
			select name into : ordvar separated by ' ' from vars where type=1;
		quit;
  %end;
  data vars; set vars (where=(type=1)); run;
  proc sql noprint;
    	select count(*) into : ck from vars;
  quit;
  %put number of ordinal variables = &ck;

%** If there are ordinal variables to process **;
  %if &ck > 0 %then %do;
    %** Get missing counts to check variable validity **;
    proc means data=&insdn noprint ;
      var &ordvar;
      output out=nmiss (drop=_type_ _freq_) nmiss= ;
    run;
    proc transpose data=nmiss out=nmiss  name=name; run;
    %** Combine **;
    proc sql;
      create table vars2 as
	  select a.*, b.col1 as miss_cnt
	  from vars a
      left join nmiss b on a.name = b.name;
    quit;
	data vars2; set vars2;
	  length new_var $32.;
	  length sign $3.;
	  length relationship $6.;
      length max_cnt unique_cnt 8.;
    run;

	%if %upcase(&transformationO) = Y %then %let maxlen=29;
	  %else %let maxlen=32;
    
	%** Create macro variables to control processing **;
    data _null_; set vars2 (where=(length(strip(name)) + length(strip("&prefix")) <= &maxlen)) end=eof;
      if eof then call symput("VARCNT",_N_);
      call symput("IV"|| trim(left(put(_N_,4.)))  ,name);
	  call symput("miss"|| trim(left(put(_N_,4.)))  ,miss_cnt);
    run;
	%** Loop through for each variable **;
    %do _I_ = 1 %to &varcnt;
	  %put processing ordinal variable &&iv&_I_;
      %if %sysevalf((&nobs*&missrate)<&&miss&_I_) %then %put missing count is &&miss&_I_ which is too high;
        %else %do;
          %put missing count is &&miss&_I_ ;
		  %pord(&insdn,&&iv&_I_,&&miss&_I_);
		%end;
    %end;
	%** Create permanent variable summary file **;
	data out.CE2_ordinal_vars; set vars2; run;
	%** Write out list of variables to keep to code file **;
	%let rck = 0;
	data _NULL_; set out.CE2_ordinal_vars (where=(missing(new_var)=0)) end=eof;
	  m = mod(_N_,7);
      length v $256.;
      retain v ;
      if m = 1 then v = new_var;
        else v = strip(v) || " " || strip(new_var);
	  FILE  "&Path_output.CE2_Ordinal_Var_Recode.txt" mod;
	  if _N_=1 then do;
        PUT ' ' '0d'x;
		  PUT ' ' '0d'x;
		PUT '%LET KEEPLIST_O =';
      end;
      if m = 0 or eof then PUT v '0d'x;
	  if eof then PUT ';';
	  if eof then call symputx("rck" ,_N_);
    run;
	%put number of recoded ordinal variables = &rck;

	%* Capture old and new variable name;
	data tmp (keep=orig_var orig_label variable);
	  length variable orig_var $32.;
	  length orig_label $256.;
      set out.CE2_ordinal_vars (where=(missing(variable)=0) rename=(name=orig_var label=orig_label new_var=variable));
	run;
	data keep_vars; set keep_vars tmp; run;

	proc datasets nolist; delete vars vars2 tmp; quit; run;
  %end;
%** If there are no ordinal variables to process **;
  %else %do;
    data _null_;
      FILE  "&Path_output.CE2_Ordinal_Var_Recode.txt" mod;
      PUT ' ' '0d'x;
      PUT '%LET KEEPLIST_O = ;' '0d'x;
    run;
  %end;
%mend ord_cntl;

** Continuous variable processing macro **;
%macro pcont(insdn,var,nmiss);
  %* check if dep_var has the constant value in nonmissing part;
  proc sql noprint;
    select min(&dep_var)=max(&dep_var) into : cons_dep from &insdn
    where missing(&var)=0;
  quit;

  %if &cons_dep = 1 %then %put &var has constant dependent value in nonmissing part and will be excluded;
  %else %do;
    %pnum(&insdn,&var,&nmiss,C);
	%** Profiling **;
    %if %upcase(&profiling) = Y %then %do;
      %prof2(&insdn,&var);
	  %prof3(C);
    %end;
  %end;

%mend pcont;
 
** Continuous processing control macro **;
%macro cont_cntl(insdn);
%** Initialize code file **;
  data _null_;
    FILE  "&Path_output.CE2_Continuous_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                      CE2_CONTINUOUS_VAR_RECODE.TXT                      ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
  run;
%** Get variables to process **;
  proc contents data=&insdn (keep=&contvar)
      out=vars (keep=varnum name type length label)
      noprint;
  run;
  proc sql noprint;
		select max(type) into : ck_type from vars;
  quit;
  %if &ck_type = 2 %then %do;
    	proc sql noprint; select name into : bad separated by ',' from vars where type = 2; quit;
		%put The following character variables were dropped from the continuous variable list:;
		%put &bad;
		proc sql noprint; 
			select name into : contvar separated by ' ' from vars where type=1;
		quit;
  %end;
  data vars; set vars (where=(type=1)); run;
  proc sql noprint;
    select count(*) into : ck from vars;
  quit;
  %put number of continuous variables = &ck;

%** If there are continuous variables to process **;
%** If there are continuous variables to process **;
  %if &ck > 0 %then %do;
    %** Get statistics to check variable validity **;
    proc means data=&insdn noprint ;
      var &contvar;
      output out=nmiss nmiss= ;
    run;
    proc transpose data=nmiss out=nmiss name=name; run;
    %** Combine **;
    proc sql;
      create table vars2 as
	  select a.*, b.col1 as miss_cnt
	  from vars a
      left join nmiss b on a.name = b.name;
    quit;
	%** Add new variables to summary file **;
	data vars2; set vars2;
	  length P&p_lo P&p_hi 8.;
	  length new_var $32.;
	  length sign $3.;
	  length relationship $6.;
	run;
	proc datasets nolist; delete nmiss; quit; run;

	%if %upcase(&transformationC) = Y %then %let maxlen=29;
	  %else %let maxlen=32;
    
	%** Create macro variables to control processing **;
    data _null_; set vars2 (where=(length(strip(name)) + length(strip("&prefix")) <= &maxlen)) end=eof;
      if eof then call symputx("VARCNT",_N_);
      call symputx("IV"|| trim(left(put(_N_,4.)))  ,name);
	  call symputx("miss"|| trim(left(put(_N_,4.)))  ,miss_cnt);
    run;
	%** Loop through for each variable **;
    %do _I_ = 1 %to &varcnt;
	  %put processing continuous variable &&iv&_I_;
      %if %sysevalf((&nobs*&missrate)<&&miss&_I_) %then %put missing count is &&miss&_I_ which is too high;
        %else %do;
          %put missing count is &&miss&_I_ ;
			 proc means data=&insdn noprint ;
      		var &&iv&_I_;
      		output out=st P&p_lo=P_lo P&p_hi=P_hi ;
    		 run;
			 data _null_; set st;
	  			call symputx("lo&_I_" ,P_lo);
	  			call symputx("hi&_I_",P_hi);
    		run;
			 data vars2; set vars2;
	  			if name = "&&iv&_I_" then do;
        			P&p_lo = &&lo&_I_;
        			P&p_hi = &&hi&_I_;
		      end; 
			 run;
		  %if %sysevalf(&&lo&_I_ = &&hi&_I_) %then %put variable has constant value and will be excluded;
		    %else %do;
			  %pcont(&insdn,&&iv&_I_,&&miss&_I_);
			%end;
		%end;
    %end;

	%** Create permanent variable summary file **;
	data out.CE2_continuous_vars; set vars2; run;
	%** Write out list of variables to keep to code file **;
	%let rck = 0;
	data _NULL_; set out.CE2_continuous_vars (where=(missing(new_var)=0)) end=eof;
	  m = mod(_N_,7);
      length v $256.;
      retain v ;
      if m = 1 then v = new_var;
        else v = strip(v) || " " || strip(new_var);
	  FILE  "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
	  if _N_=1 then do;
        PUT ' ' '0d'x;
		  PUT ' ' '0d'x;
		PUT '%LET KEEPLIST_C =';
      end;
      if m = 0 or eof then PUT v '0d'x;
	  if eof then PUT ';';
	  if eof then call symputx("rck" ,_N_);
    run;
	 %if &rck = 0 %then %do;
    	data _null_;
      	FILE  "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
      	PUT ' ' '0d'x;
      	PUT '%LET KEEPLIST_C = ;' '0d'x;
    	run;
  	 %end;
	%put number of recoded continuous variables = &rck;

	%* Capture old and new variable name;
	data tmp (keep=orig_var orig_label variable);
	  length variable orig_var $32.;
	  length orig_label $256.;
      set out.CE2_continuous_vars (where=(missing(variable)=0) rename=(name=orig_var label=orig_label new_var=variable));
	run;
	data keep_vars; set keep_vars tmp; run;

	proc datasets nolist; delete vars vars2 tmp; quit; run;
  %end;
%** If there are no continuous variables to process **;
  %else %do;
    data _null_;
      FILE  "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
      PUT ' ' '0d'x;
      PUT '%LET KEEPLIST_C = ;' '0d'x;
    run;
  %end;
%mend cont_cntl;

** Master EDA and Recode macro **;
%macro CE_EDA_Recode(insdn=out.CE1_Resampled);
%** Extract working portion of file **;
data workfile; set &insdn (where=(mod_val_test ^= 3)); run;

%* Get global numbers;
%global nobs overall_avg minbinn;
proc sql noprint;
  select count(*) into : nobs from workfile;
  select mean(&dep_var) into : overall_avg from workfile;
quit;
data _null_;
  call symputx("minbinn",max(&minbinnc,(&nobs*&minbinnp)));     ** minimum bin size for nominal;
run;

%** Initialize profiling dataset **;
%if %upcase(&profiling) = Y %then %do;
  data out.CE2_profile; set _NULL_;
    length type $10.;
    length variable $32.;
	length label $256.;
    length category $256.;
	format count comma8.;
	format percent percent8.2;
	%if %upcase(&binary_dv) = Y %then %do;
      format Average_DV percent8.2;
	%end;
	%else %do;
	  format Average_DV 12.2;
	%end;
	format index 8.0;
	length star $8.;
  run;
%end;

%** Initialize variable name lookup dataset **;
data keep_vars; set _NULL_;
  length orig_var $32.;
  length orig_label $256.;
  length variable $32.;
run;

%if %symexist(binvar) %then %do;
  %bin_cntl(workfile);
%end;
%else %do;
  data _null_;
    FILE  "&Path_output.CE2_Binary_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                        CE2_BINARY_VAR_RECODE.TXT                        ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
	PUT ' ' '0d'x;
    PUT '%LET KEEPLIST_B = ;' '0d'x;
  run;
%end;

%if %symexist(nomvar) %then %do;
  %nom_cntl(workfile);
%end;
%else %do;
  data _null_;
    FILE  "&Path_output.CE2_Nominal_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                        CE2_NOMINAL_VAR_RECODE.TXT                       ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
    PUT ' ' '0d'x;
    PUT '%LET KEEPLIST_N = ;' '0d'x;
  run;
%end;

%if %symexist(ordvar) %then %do;
  %ord_cntl(workfile);
%end;
%else %do;
  data _null_;
    FILE  "&Path_output.CE2_Ordinal_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                       CE2_ORDINAL_VAR_RECODE.TXT                        ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
	PUT ' ' '0d'x;
    PUT '%LET KEEPLIST_O = ;' '0d'x;
  run;
%end;

%if %symexist(contvar) %then %do;
  %cont_cntl(workfile);
%end;
%else %do;
  data _null_;
    FILE  "&Path_output.CE2_Continuous_Var_Recode.txt" LRECL=256;
    PUT "*********************************************************************************;" '0d'x/
        "****                      CE2_CONTINUOUS_VAR_RECODE.TXT                      ****;" '0d'x/
	    "*********************************************************************************;" '0d'x;
	PUT ' ' '0d'x;
    PUT '%LET KEEPLIST_C = ;' '0d'x;
  run;
%end;

%** Create dataset with recoded variables **;
data out.CE2_Recoded; set &insdn;
  %inc "&Path_output.CE2_Binary_Var_Recode.txt";
  %inc "&Path_output.CE2_Nominal_Var_Recode.txt";
  %inc "&Path_output.CE2_Ordinal_Var_Recode.txt";
  %inc "&Path_output.CE2_Continuous_Var_Recode.txt";
  keep &keep_list &KEEPLIST_B &KEEPLIST_N &KEEPLIST_O &KEEPLIST_C;
run;

%** Finalize variable name lookup dataset **;
proc contents data=out.CE2_Recoded out=vars (keep=name label) noprint; run;
proc sql;
  create table out.CE2_vars as
  select a.variable, b.label, a.orig_var, a.orig_label
  from keep_vars a, vars b
  where a.variable = b.name
  order by a.variable;
quit;

%** Create EDA report in Excel;
ods listing close;
ods Tagsets.ExcelxP body="&Path_output.CE2_EDA_report.xls" style=sasweb;

%if %sysfunc(exist(out.CE2_binary_vars)) %then %do;
	proc sql noprint;
		select max(length(name)) into: l1 from out.CE2_binary_vars;
		select max(length(label)) into: l2 from out.CE2_binary_vars;
		select max(length(new_var)) into: l3 from out.CE2_binary_vars;
	quit;
	ods tagsets.excelxp options(sheet_name="Binary Variables" absolute_column_width="&l1,8,8,&l2,&l3,8,8" 
		Frozen_Headers='Yes' Frozen_RowHeaders='1');
	proc print data=out.CE2_binary_vars noobs;
		var name type length label new_var;
		var miss_cnt ones_cnt / style={tagattr='format:#,##0'};
	run;
%end;
%if %sysfunc(exist(out.CE2_nominal_vars)) %then %do;
	proc sql noprint;
		select max(length(name)) into: l1 from out.CE2_nominal_vars;
		select max(length(label)) into: l2 from out.CE2_nominal_vars;
		select max(length(new_vars)) into: l3 from out.CE2_nominal_vars;
	quit;
	ods tagsets.excelxp options(sheet_name="Nominal Variables" absolute_column_width="&l1,8,8,&l2,&l3,8,8,8"
		Frozen_Headers='Yes' Frozen_RowHeaders='1');
	proc print data=out.CE2_nominal_vars noobs;
		var name type length label new_vars;
		var miss_cnt unique_cnt max_cnt / style={tagattr='format:#,##0'};
	run;
%end;
%if %sysfunc(exist(out.CE2_ordinal_vars)) %then %do; 
	proc sql noprint;
		select max(length(name)) into: l1 from out.CE2_ordinal_vars;
		select max(length(label)) into: l2 from out.CE2_ordinal_vars;
		select max(length(new_var)) into: l3 from out.CE2_ordinal_vars;
	quit;
	ods tagsets.excelxp options(sheet_name="Ordinal Variables" 
		%if %upcase(&stdmethodC) ^= NO %then %do; 
			%if %upcase(&binary_dv)=Y %then %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,8,8,9,8,8,8,8" %end;
			%else %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,8,8,9,8,8,8" %end;
		%end;
		%else %do; 
			%if %upcase(&binary_dv)=Y %then %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,9,8,8,8,8" %end;
			%else %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,9,8,8,8" %end;
		%end;
		Frozen_Headers='Yes' Frozen_RowHeaders='1');
	proc print data=out.CE2_ordinal_vars noobs;
		var name type length label new_var relationship;
		var miss_cnt unique_cnt max_cnt / style={tagattr='format:#,##0'};
		var var_lb var_ub var_median / style={tagattr='format:0.00'};
		%if %upcase(&stdmethodO) ^= NO %then %do; var Location Scale / style={tagattr='format:0.00'}; %end;
		var miss_impute / style={tagattr='format:0.0000'};
		%if %upcase(&Binary_dv)= Y %then %do;
			var Concordant / style={tagattr='format:0.00'};
			var CValue Chisq PValue / style={tagattr='format:0.0000'};
		%end;
		%else %do;
			var RSquare FValue PValue / style={tagattr='format:0.0000'};
		%end;
	run;
%end;
%if %sysfunc(exist(out.CE2_continuous_vars)) %then %do; 
	proc sql noprint;
		select max(length(name)) into: l1 from out.CE2_continuous_vars;
		select max(length(label)) into: l2 from out.CE2_continuous_vars;
		select max(length(new_var)) into: l3 from out.CE2_continuous_vars;
	quit;
	ods tagsets.excelxp options(sheet_name="Continuous Variables" 
		%if %upcase(&stdmethodC) ^= NO %then %do; 
			%if %upcase(&binary_dv)=Y %then %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,8,8,9,8,8,8,8" %end;
			%else %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,8,8,9,8,8,8" %end;
		%end;
		%else %do; 
			%if %upcase(&binary_dv)=Y %then %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,9,8,8,8,8" %end;
			%else %do;
				absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,9,8,8,8" %end;
		%end;
		Frozen_Headers='Yes' Frozen_RowHeaders='1');
	proc print data=out.CE2_continuous_vars noobs;
		var name type length label new_var relationship;
		var miss_cnt / style={tagattr='format:#,##0'};
		var P&p_lo P&p_hi var_lb var_ub var_median / style={tagattr='format:0.00'};
		%if %upcase(&stdmethodC) ^= NO %then %do; var Location Scale / style={tagattr='format:0.00'}; %end;
		var miss_impute / style={tagattr='format:0.0000'};
		%if %upcase(&Binary_dv)= Y %then %do;
			var Concordant / style={tagattr='format:0.00'};
			var CValue Chisq PValue / style={tagattr='format:0.0000'};
		%end;
		%else %do;
			var RSquare FValue PValue / style={tagattr='format:0.0000'};
		%end;
	run;
%end;

ods Tagsets.ExcelxP close;
ods listing;

proc corr data=out.CE2_Recoded noprint outp=corr;
	var &dep_var;
	with &PREFIX.:;
run;
proc sql;
	create table out.CE2_corr as
	select a.name as variable, a.label, b.&dep_var as correlation
	from vars a, corr b
	where a.name = b._name_
	order by abs(b.&dep_var) desc;
quit;

proc sql noprint;
	select max(length(variable)) into: l1 from out.CE2_corr;
	select max(length(label)) into: l2 from out.CE2_corr;
quit;

ods listing close;
ods Tagsets.ExcelxP body="&Path_output.CE2_Correlation_report.xls" style=sasweb;
ods tagsets.excelxp options(sheet_name="Correlations" absolute_column_width="&l1,&l2,8" Frozen_Headers='Yes' Frozen_RowHeaders='1');
proc print data=out.CE2_corr noobs;
	var variable label;
	var correlation / style={tagattr='format:0.0000'};
run;
ods Tagsets.ExcelxP close;
ods listing;

%if %upcase(&profiling) = Y %then %do;
proc sql noprint;
	select max(length(variable)) into: l1 from out.CE2_profile;
	select max(length(label)) into: l2 from out.CE2_profile;
	select max(length(category)) into: l3 from out.CE2_profile;
quit;

%** Create Profile report in Excel;
%* Order file;
proc sql;
	create table ord as
	select variable, max(abs(index))-min(abs(index)) as diff
	from out.CE2_profile
	group by variable
	order by diff desc;
quit;
data ord; set ord;
	sord = _n_;
run;
data prof; set out.CE2_profile; s2 = _n_; run;
proc sql;
	create table prof2 as
	select a.*
	from prof a, ord b
	where a.variable = b.variable
	order by b.sord, a.s2;
quit;
data out.CE2_profile; set prof2 (drop= s2); run;
proc datasets nolist; delete prof prof2; quit; run;

ods listing close;
ods Tagsets.ExcelxP body="&Path_output.CE2_Profile_report.xls" style=sasweb;

ods tagsets.excelxp options(sheet_name="Binary Variables" absolute_column_width="&l1,&l2,&l3,8,8,9,8,6"
	Frozen_Headers='Yes' Frozen_RowHeaders='1');
proc print data=out.CE2_profile noobs; 
	where type = "Binary";
	var variable label category;
	var count / style={tagattr='format:#,##0'};
	var percent Average_DV index star;
run;
ods tagsets.excelxp options(sheet_name="Nominal Variables" absolute_column_width="&l1,&l2,&l3,8,8,9,8,6"
	Frozen_Headers='Yes' Frozen_RowHeaders='1');
proc print data=out.CE2_profile noobs; 
	where type = "Nominal";
	var variable label category;
	var count / style={tagattr='format:#,##0'};
	var percent Average_DV index star;
run;
ods tagsets.excelxp options(sheet_name="Ordinal Variables" absolute_column_width="&l1,&l2,&l3,8,8,9,8,6"
	Frozen_Headers='Yes' Frozen_RowHeaders='1');
proc print data=out.CE2_profile noobs; 
	where type = "Ordinal";
	var variable label category;
	var count / style={tagattr='format:#,##0'};
	var percent Average_DV index star;
run;
ods tagsets.excelxp options(sheet_name="Continuous Variables" absolute_column_width="&l1,&l2,&l3,8,8,9,8,6"
	Frozen_Headers='Yes' Frozen_RowHeaders='1');
proc print data=out.CE2_profile noobs; 
	where type = "Continuous";
	var variable label category;
	var count / style={tagattr='format:#,##0'};
	var percent Average_DV index star;
run;

ods Tagsets.ExcelxP close;
ods listing;
%end;

proc datasets nolist; delete workfile corr vars keep_vars; quit; run;
%MEND;              

**********************************************************************************************************************;
*** 3.Variable reduction and ranking ***;
**********************************************************************************************************************;
** Master variable reduction macro **;
%macro CE_Var_Redu(insdn=out.CE2_Recoded);
	%* Check if the user has specified a fast run;
	%if %symexist(fast_opt) %then %do;
		%if %upcase(&fast_opt) = Y %then %do;
			%let sources = 1;
			%let univ_reg = N;
			%let correlation = N;
			%let principal = N;
			%let cluster = N;
			%if %upcase(&binary_dv) = Y %then %do;
				%let regression = N;
				%let logistic = Y;
			%end;
			%else %do;
				%let regression = Y;
				%let logistic = N;
			%end;
			%let information = N;
			%let ind_correlation = N;
			%let ind_dv_corr=N;
		%end;
	%end;

  %** Create working file;
  proc sql noprint;
    select count(*) into : nobs from &insdn where mod_val_test^=3;
  quit;
  %put observation count is &nobs;
  data workfile; set &insdn (where=(mod_val_test^=3));
    %if &nobs > &samplesize %then %do; if uniform(131071)<=1.0*&samplesize/&nobs; %end; 
  run;

  %** Get independent variables;
  proc contents data=workfile (drop=&keep_list)
    out=vars (keep=name label) noprint;
  run;
  proc sql noprint;
    select count(*) into : varnum from vars;
  quit;
  %put Number of variables is &varnum;

  %** Reorder variables to avoid bias when working in groups;
  data vars; set vars; rannum = ranuni(274923); run;
  proc sort data = vars; by rannum; run;

  %** Using Univariate Regression method **;
  %if %upcase(&univ_reg)=Y %then %do;
    %do i=1 %to &varnum;
      data _null_; set vars;
	    if _n_=&i;
	    call symputx('curvar',name);
      run;
      %Univ_Reg(workfile,&curvar);
	  %if &i = 1 %then %do;
        data outsdnuniv; length variable $32.; set univ_tmp; run;
	  %end;
	  %else %do;
	    data outsdnuniv; set outsdnuniv univ_tmp; run;
	  %end;
    %end;
	proc datasets library=work nolist; delete univ_tmp; quit; run;
    %* Finalize file;
    data outsdnuniv; set outsdnuniv;
      if PValue <= &maxpuni then univ_flag = 1; else univ_flag = 0;
    run;
    proc sort data = outsdnuniv; by variable; run;
  %end;

  %** Using correlation method **;
  %if %upcase(&correlation)=Y %then %do;
    proc corr data=workfile noprint outp=outsdncorr;
      var &dep_var;
      with &PREFIX.:;
	  %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end; 
    run;
    %* Finalize file;
    data outsdncorr (drop=_type_);
	  length variable $32.;
      set outsdncorr (where=(_type_='CORR') rename=(_name_=variable &dep_var=corr));
      if abs(corr) >= &corrcut then corr_flag=1; else corr_flag=0;
    run;
   	proc sort data= outsdncorr; by variable; run;
  %end;

  %** Using factor/principal analysis method **;
  %if %upcase(&principal)=Y %then %do;
    %if &varnum <= &nprin %then %do; 
      %let nprin = %eval(&varnum-1); 
    %end;
	%** Split variables into manageable groups;
	%let group=%sysfunc(max(%eval(&varnum/(&nprin*5)),1));
    data tmp; set vars; rank=ceil(_n_/(&varnum/&group)); run;
	%** Process variables by group;
	%do i = 1 %to &group;
	  proc sql noprint; select name into : regvlist separated by ' ' 
        from tmp where rank = &i; 
      quit;
      %single_group_prin(&insdn,prin_tmp,&regvlist);
	  %if &i = 1 %then %do;
        data outsdnprin; length variable $32.; set prin_tmp; run;
	  %end;
	  %else %do;
	    data outsdnprin; set outsdnprin prin_tmp; run;
	  %end;
    %end;
	proc datasets library=work nolist; delete prin_tmp tmp; quit; run;
    %* Finalize file;
    %if &group > 1 %then %do;
      proc sql noprint; select variable into : regvlist separated by ' ' 
        from outsdnprin where Factor >= &minprin; 
      quit;
	  %single_group_prin(&insdn,outsdnprin,&regvlist);
    %end;
    data outsdnprin; set outsdnprin;
      if Factor >= &minprin then prin_flag=1; else prin_flag=0;
    run;
	proc sort data = outsdnprin; by variable; run;
  %end;

  %** Using Clustering method **;
  %if %upcase(&cluster)=Y %then %do;
    %if &varnum <= &maxc %then %do; 
      %let maxc = %eval(&varnum-1); 
    %end;
	%** Split variables into manageable groups;
	%let group=%sysfunc(max(%eval(&varnum/(&maxc*5)),1));
    data tmp; set vars; rank=ceil(_n_/(&varnum/&group)); run;
	%** Process variables by group;
	%do i = 1 %to &group;
	  proc sql noprint; select name into : regvlist separated by ' ' 
        from tmp where rank = &i; 
      quit;
      %single_group_clus(&insdn,clus_tmp,&regvlist);
	  %if &i = 1 %then %do;
        data outsdnclus; length variable $32.; set clus_tmp; run;
	  %end;
	  %else %do;
	    data outsdnclus; set outsdnclus clus_tmp; run;
	  %end;
    %end;
	proc datasets library=work nolist; delete clus_tmp tmp; quit; run;
    %* Finalize file;
    %if &group > 1 %then %do;
      proc sql noprint; select variable into : regvlist separated by ' ' 
        from outsdnclus where RSquareRatio <= &maxratio; 
      quit;
	  %single_group_clus(&insdn,outsdnclus,&regvlist);
    %end;
    data outsdnclus; set outsdnclus;
      if RSquareRatio <= &maxratio then clus_flag=1; else clus_flag=0;
    run;
    proc sort data = outsdnclus; by variable; run;
  %end;

  %** Using regression method **;
  %if %upcase(&regression)=Y   %then %do;
	%** Split variables into manageable groups;
	%let group=%sysfunc(ceil(&varnum/100));
    data tmp; set vars; rank=ceil(_n_/(&varnum/&group)); run;
	%** Process variables by group;
	%do i = 1 %to &group;
	  proc sql noprint; select name into : regvlist separated by ' ' 
        from tmp where rank = &i; 
      quit;
      %single_group_reg(&insdn,reg_tmp,&regvlist);
	  %if &i = 1 %then %do;
        data outsdnreg; length variable $32.; set reg_tmp; run;
	  %end;
	  %else %do;
	    data outsdnreg; set outsdnreg reg_tmp; run;
	  %end;
    %end;
	proc datasets library=work nolist; delete reg_tmp tmp; quit; run;
    %* Finalize file;
    %if &group > 1 %then %do;
      proc sql noprint; select variable into : regvlist separated by ' ' 
        from outsdnreg; 
      quit;
	  %single_group_reg(&insdn,outsdnreg,&regvlist);
    %end;
    proc sort data = outsdnreg; by variable; run;
  %end;

  %** Using logistic method **;
  %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y  %then %do;
    %** Split variables into manageable groups;
	%let group=%sysfunc(ceil(&varnum/100));
    data tmp; set vars; rank=ceil(_n_/(&varnum/&group)); run;
	%** Process variables by group;
	%do i = 1 %to &group;
	  proc sql noprint; select name into : regvlist separated by ' ' 
        from tmp where rank = &i; 
      quit;
      %single_group_log(&insdn,log_tmp,&regvlist);
	  %if &i = 1 %then %do;
        data outsdnlog; length variable $32.; set log_tmp; run;
	  %end;
	  %else %do;
	    data outsdnlog; set outsdnlog log_tmp; run;
	  %end;
    %end;
	proc datasets library=work nolist; delete log_tmp tmp; quit; run;
    %* Finalize file;
    %if &group > 1 %then %do;
      proc sql noprint; select variable into : regvlist separated by ' ' 
        from outsdnlog; 
      quit;
	  %single_group_log(&insdn,outsdnlog,&regvlist);
    %end;
    proc sort data = outsdnlog; by variable; run;
  %end;

  %** Using information value method **;
  %if %upcase(&information)=Y %then %do;
    %do i=1 %to &varnum;
      data _null_; set vars;
	    if _n_=&i;
	    call symputx('curvar',name);
      run;
      %Info_Val_Var(workfile,&curvar,&varnum);
	  %if &i = 1 %then %do;
        data outsdninfv; length variable $32.; set infv_tmp; run;
	  %end;
	  %else %do;
	    data outsdninfv; set outsdninfv infv_tmp; run;
	  %end;
    %end;
	proc datasets library=work nolist; delete infv_tmp; quit; run;
    %* Finalize file;
    data outsdninfv; set outsdninfv;
      if infv >= &infvcut then infv_flag = 1; else infv_flag = 0;
    run;
    proc sort data = outsdninfv; by variable; run;
  %end;
  
  %* Get basic metrics for variables;
  proc means data=&insdn noprint ;
    var &prefix: ;
    output out=cnt n= ;
    output out=min min= ;
    output out=max max= ;
	output out=mean mean= ;
  run;
  data stats;
    set cnt (in=a) min (in=b) max (in=c) mean (in=d);
	length var $8;
	if a then var = "Count";
	  else if b then var = "Minimum";
	  else if c then var = "Maximum";
	  else if d then var = "Mean";
    drop _freq_ _type_;
  run;
  proc transpose data=stats out=stats name=Variable label=Label; id var; run;
  proc sort data=stats; by variable; run;
  data stats; length variable $32.; set stats; run;

  %*combine results from all methods together;
  data out.CE3_Var_Redu;
    merge stats
      %if %upcase(&univ_reg)=Y %then %do; outsdnuniv(in=tuniv) %end;
      %if %upcase(&correlation)=Y %then %do; outsdncorr(in=tcorr) %end;
      %if %upcase(&principal)=Y %then %do; outsdnprin(in=tprin) %end;
      %if %upcase(&cluster)=Y %then %do; outsdnclus(in=tclus) %end;
      %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; outsdnlog(in=tlog) %end;
      %if %upcase(&regression)=Y %then %do; outsdnreg(in=treg) %end;
      %if %upcase(&information)=Y %then %do; outsdninfv(in=tinfv) %end;
    ;
    by variable;
      %if %upcase(&univ_reg)=Y %then %do; if tuniv and univ_flag=1 then univsource=1; else univsource=0; drop univ_flag; %end;
      %if %upcase(&correlation)=Y %then %do; if tcorr and corr_flag=1 then corrsource=1; else corrsource=0; drop corr_flag; %end;
      %if %upcase(&principal)=Y %then %do; if tprin and prin_flag=1 then prinsource=1; else prinsource=0; drop prin_flag; %end;
      %if %upcase(&cluster)=Y %then %do; if tclus and clus_flag=1then clussource=1; else clussource=0; drop clus_flag; %end;
      %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; if tlog then logsource=1; else logsource=0; %end;
      %if %upcase(&regression)=Y %then %do; if treg then regsource=1; else regsource=0; %end;
      %if %upcase(&information)=Y %then %do; if tinfv and infv_flag=1 then infvsource=1; else infvsource=0; drop infv_flag; %end;

    num_sources=sum(of 
      %if %upcase(&univ_reg)=Y %then %do; univsource  %end;
      %if %upcase(&correlation)=Y %then %do; corrsource  %end;
      %if %upcase(&principal)=Y %then %do; prinsource  %end;
      %if %upcase(&cluster)=Y %then %do; clussource  %end;
      %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; logsource %end;
      %if %upcase(&regression)=Y %then %do; regsource %end;
      %if %upcase(&information)=Y %then %do; infvsource %end;
    );
  run;
  proc sort data=out.CE3_Var_Redu; by descending num_sources %if %upcase(&information)=Y %then %do; descending infv %end; ;run;

  %** Exclude highly correlated variables **;
  %if %upcase(&ind_correlation)=Y %then %do;
    %Maxx_Corr(workfile);
  %end;

  %** Exclude variables highly correlated to dependent variable **;
  %if %upcase(&ind_dv_corr)=Y %then %do;
    proc sql;
	   create table tmp as
		select a.*, b.correlation as dv_corr,
		  case when b.correlation > &max_dv_corr then 'Y' end as drop_dv_corr
		from out.CE3_Var_Redu a
		left join out.CE2_corr b on a.variable = b.variable;
	 quit;
	 data out.CE3_Var_Redu; set tmp; run;
	 proc datasets library=work nolist; delete tmp; quit; run;
  %end;

  %** Create report in Excel;
  ods listing close;
  ods Tagsets.ExcelxP body="&Path_output.CE3_Var_Redu Results.xls" style=sasweb;
  ods tagsets.excelxp options(sheet_name="Variables");
    proc print data=out.CE3_Var_Redu (drop=
	           %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; LogStep %end;
			   %if %upcase(&regression)=Y %then %do; RegStep %end;
			   %if %upcase(&cluster)=Y %then %do; Cluster  %end;
         ) noobs; 
    run;
  ods Tagsets.ExcelxP close;
  ods listing;

  %*Output list of top variables;
  data selected; set out.CE3_Var_Redu;
    %if %upcase(&ind_correlation)=Y %then %do; if missing(drop_corr); %end;
	 %if %upcase(&ind_dv_corr)=Y %then %do; if missing(drop_dv_corr); %end;
    if num_sources>=&sources 
      %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; or logsource=1 %end;
      %if %upcase(&regression)=Y %then %do; or regsource=1 %end; ;
  run;
  %let rck = 0;
  data _null_; set selected end=eof;
    m = mod(_N_,7);
    length v $256.;
    retain v ;
    if m = 1 then v = variable;
      else v = strip(v) || " " || strip(variable);
	FILE  "&Path_output.CE3_Varlist_redu.txt" lrecl=256;
	if _N_=1 then PUT '%let varlist_redu =';
    if m = 0 or eof then PUT v '0d'x;
	if eof then PUT ';';
	if eof then call symputx("rck" ,_N_);
  run;
  %put number of selected variables = &rck;
  proc datasets nolist; delete stats cnt min max mean selected workfile vars
      %if %upcase(&univ_reg)=Y %then %do; outsdnuniv %end;
      %if %upcase(&correlation)=Y %then %do; outsdncorr %end;
      %if %upcase(&principal)=Y %then %do; outsdnprin %end;
      %if %upcase(&cluster)=Y %then %do; outsdnclus %end;
      %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; outsdnlog %end;
      %if %upcase(&regression)=Y %then %do; outsdnreg %end;
      %if %upcase(&information)=Y %then %do; outsdninfv %end;
  ; quit; run;

%mend;

** Univ_Reg: using univariate regression to do variable reduction **;
%macro Univ_Reg(insdn,var);
  data univ_tmp; set _NULL_; run;

  %PUT ***Univ_Reg STEP, CURRENT VARIABLE: &var***;
  %if %upcase(&Binary_dv) = Y %then %do; 
    %* Run univariate logistic regression;
    ods listing close;                                                                             
    ods output parameterestimates=parm association=fitstat1;                                          
	proc logistic data=&insdn desc namelen=32;
 	  model &dep_var=&var;
	  %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end; 
    run;
  	ods listing; 
	%if %sysfunc(exist(fitstat1)) %then %do; 
	  data fitstat1 (keep=CC_RSQ); set fitstat1(keep=cvalue1 obs=1);
	    CC_RSQ = input(cvalue1,best8.);
	  run;
	%end;
	%else %do;
	  data fitstat1; CC_RSQ = .; output; run;
	%end;
    data parm (keep=ProbChiSq sign rename=(ProbChiSq=PValue)); 
	  length sign $6. ;
      set parm (firstobs=2 obs=2 keep=variable Estimate ProbChiSq);
      if Estimate>0 then sign='(+)'; else sign='(-)';
  	run;
  %end;
  %else %do;
	%* Run univariate simple regression;
	ods listing close;
    ods output SelParmEst=parm SelectionSummary=fitstat1;
	proc reg data=&insdn;
	  model &dep_var=&var /selection=forward MAXSTEP=1 slentry=0.999;
	  %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end;
	run;
 	ods listing;
	%if %sysfunc(exist(fitstat1)) %then %do; 
	  data fitstat1 (rename=(ModelRsquare=CC_RSQ)); set fitstat1(keep=ModelRsquare obs=1); run;
	%end;
	%else %do;
	  data fitstat1; CC_RSQ = .; output; run;
	%end;
	data parm (keep=ProbF sign rename=(ProbF=PValue)); 
	  length sign $6. ;
      set parm (firstobs=2 obs=2 keep=variable Estimate ProbF);
      if Estimate>0 then sign='(+)'; else sign='(-)';
  	run;
  %end;
  %* Combine;
  data univ_tmp;
	length variable $32.;
    merge parm fitstat1;
    variable="&var";
  run;
  %* Clean up;
  proc datasets library=work nolist; delete fitstat1 parm; quit; run;
%mend;

** Single_Group_Prin: using principal components to do variable reduction **;
%macro Single_Group_Prin(insdn,outsdn,vlist);
  %local j;
  data &outsdn; set _NULL_; run;
  %* Run factor / principal components analysis;
  proc factor data=&insdn out=tmpdataa2 method=prin priors=one nfact=&nprin noprint;
    var &vlist;
    %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end;
  run;
  %* Get correlations;
  proc corr data=tmpdataa2 out=tmpdataa3 noprint;
    var Factor: ;
    with &vlist;
    %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end;
  run;
  %* Rearrange dataset;
  data tmpdataa4 (drop=_type_ _name_); set tmpdataa3 (where=(_type_='CORR')); 
    length variable $32;  
    variable=_name_;
    %do j = 1 %to &nprin;
      factor=abs(factor&j); 
	  output;
	%end;
  run;
  %* Get best correlation for each variable;
  proc sort data=tmpdataa4 (keep=variable factor); by variable descending factor; run;
  data &outsdn; set tmpdataa4;
    by variable;
    if first.variable;
  run;
  proc datasets library=work nolist; delete tmpdataa: ; quit; run;
%mend;

** Single_Group_Clus: using clustering to do variable reduction **;
%macro Single_Group_Clus(insdn,outsdn,vlist);
  data &outsdn; set _NULL_; run;
  %* Run cluster procedure;
  ods listing close;
  ods output rsquare=&outsdn;
  proc varclus data=&insdn (keep=&vlist &weight) minc=&maxc maxc=&maxc short;
    var &vlist;
    %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end;
  run;
  ods listing;

  %* Clean up file;
  data &outsdn (keep=variable cluster RSquareRatio); 
    length variable $32.;
    set &outsdn; 
    retain clustertemp;
    if cluster ne ' ' then  clustertemp=cluster;
      else cluster=clustertemp;
  run;
%mend;

** Single_Group_Reg: using linear regression to do variable reduction **;
%macro Single_Group_Reg(insdn,outsdn,vlist);
  data &outsdn; set _NULL_; run;
  ods listing close;
  ods output   SelectionSummary=&outsdn;
  proc reg data=&insdn;
    model &dep_var = &vlist /selection=forward slentry=&alphareg;
    %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end;
  run;
  ods listing;

  data &outsdn (keep=variable RegStep RegPValue); 
    length variable $32.;
    set  &outsdn (rename=(VarEntered = variable step=RegStep ProbF = RegPValue));
  run;
%mend;

** Single_Group_Log: using logistic regression to do variable reduction **;
%macro Single_Group_Log(insdn,outsdn,vlist);
  data &outsdn; set _NULL_; run;
  ods listing close;
  ods output   ModelBuildingSummary=&outsdn;
  proc logistic data=&insdn desc namelen=32;
    model &dep_var = &vlist /selection=forward slentry=&alphalog;
    %if %upcase(&redu_weight) = Y %then %do; weight &weight; %end;
  run;
  ods listing;

  data &outsdn (keep=variable LogStep LogPValue); 
    length variable $32.;
    set  &outsdn (rename=(EffectEntered = variable step=LogStep ProbChiSq = LogPValue));
  run;
%mend;

** Info_Val_Var: using information value to do variable reduction **;
%macro Info_Val_Var(insdn,var,set_size);
  data infv_tmp; set _NULL_; run;
  %* Check number of unique values;
  proc sql noprint;
    select count(distinct(&var)) into : unq from &insdn;
  quit;

  %if %eval(&unq>&decile) %then %do;
    %if %upcase(&redu_weight) = Y %then %do;
	  %* Create bins if there is weighting;
      proc sql noprint;
        select sum(&weight) into : cumwgt from &insdn;
      quit;
      proc sort data=&insdn (keep=&dep_var &var &weight) out=tmp; by &var; run;
      data tmp (drop=rank); set tmp;
        retain rank;
        rank + &weight;
        bin = (floor(rank*&decile/(&cumwgt+1)));
      run;
	  %* Summarize data to bins;
	  proc summary data=tmp;
        weight &weight;
        var &var &dep_var;
        class bin;
	    output out=tmp2 (drop=_freq_) sumwgt=cnt mean(&dep_var)=mean_dv max(&dep_var)=max_dv
	                    mean(&var)=mean_var min(&var)=min_var max(&var)=max_var /noinherit;
	  run;
    %end;
    %else %do;
	  %* Create bins if there is no weighting;
      proc rank data=&insdn (keep=&dep_var &var) out=tmp groups=&decile;
        var &var;
	    ranks bin;
	  run;
	  %* Summarize data to bins;
	  proc summary data=tmp;
        var &var &dep_var;
        class bin;
	    output out=tmp2 (rename=_freq_=cnt) mean(&dep_var)=mean_dv max(&dep_var)=max_dv
	                    mean(&var)=mean_var min(&var)=min_var max(&var)=max_var /noinherit;
	  run;
	%end;
  %end;
  %else %do;
    %* Summarize to &var value when less than &decile;
    %if %upcase(&redu_weight) = Y %then %do;
      proc summary data=&insdn;
        weight &weight;
        var &dep_var;
        class &var;
	    output out=tmp2 (drop=_freq_ rename=&var=bin) sumwgt=cnt mean(&dep_var)=mean_dv max(&dep_var)=max_dv /noinherit;
	  run;
	%end;
	%else %do;
      proc summary data=&insdn;
        var &dep_var;
        class &var;
	    output out=tmp2 (rename=(_freq_=cnt &var=bin)) mean(&dep_var)=mean_dv max(&dep_var)=max_dv /noinherit;
	  run;
	%end;
	data tmp2; set tmp2;
	  mean_var = bin;
	  min_var = bin;
	  max_var = bin;
	run;
  %end;

  %* Create global metrics;
  proc sql noprint; 
    select mean_dv into : norm   from tmp2 where _type_ = 0;
    select cnt into : ntotal from tmp2 where _type_ = 0; 
	select max(max_dv) into : maxdv from tmp2 where _type_ = 0;
	select sum(cnt*mean_dv) into: totalresp from tmp2 where _type_ = 1;
  quit;

  %* Calculate information value, ks and gini statistics;
  data tmp3; set tmp2 (where=(_type_=1));
    lift_index = (mean_dv/&norm)*100; 
    retain cumresp cumtotal;
    if _n_ = 1 then do; 
      cumresp = 0; 
      cumtotal = 0; 
    end;
    cumresp = cumresp + (mean_dv*cnt);
    cumpct_resp = cumresp/&totalresp;
    cumtotal = cumtotal + cnt;
    cumpct_freq = cumtotal/&ntotal;
    cumavg = cumresp/cumtotal;
    cumindex = (cumavg/&norm)*100; 
    badrate = mean_dv/&maxdv;
    goodrate = 1 - badrate;
    badcnt = badrate * cnt;  
    goodcnt = goodrate * cnt; 
    badratio = (cnt * mean_dv)/&totalresp;
    goodratio = (cnt - badcnt)/(&ntotal-&totalresp/&maxdv); 
    drop cumresp cumtotal;
   run;  
   data infv_tmp; set tmp3 end=eof;      
	 length variable $32.;
	 variable = "&var";
     retain infv gini ks cumbad cumgood 0;
     prev_dv = lag(cumpct_resp);
     prev_pct = lag(cumpct_freq);
     if _n_ ne 1 then 
       gini = gini+2*((cumpct_freq+prev_pct)/2-(cumpct_resp+prev_dv)/2)*(cumpct_freq-prev_pct);
     if badratio+goodratio = 0 then contribution = 0;
       else if badratio*goodratio=0 then contribution=cnt*1.0/&set_size*max(badratio,goodratio);
       else contribution = (badratio-goodratio)*log(Max(badratio/goodratio,0.00001));
     infv = infv+(contribution);
     cumbad = cumbad + badratio;
     cumgood = cumgood + goodratio;
     ks = max(ks,abs(cumbad-cumgood));
     if eof then output; 
     keep variable infv gini ks;
   run;
   proc datasets library=work nolist; delete tmp: ; quit; run;
%mend;

** Maxx_Corr: eliminate variables that are highly correlated **;
%macro Maxx_Corr(insdn);
  %* Add field to variable dataset;
  data out.CE3_Var_Redu; set out.CE3_Var_Redu;
    length drop_corr $1.;
  run;

  %* Get variables to evaluate;
  data tmp; set out.CE3_Var_Redu;
    if num_sources>=&sources 
      %if %upcase(&logistic)=Y and %upcase(&Binary_dv)=Y %then %do; or logsource=1 %end;
      %if %upcase(&regression)=Y %then %do; or regsource=1 %end; ;
  run;
  proc sql noprint;
    select variable into : varlist_redu separated by ' ' from tmp;
  quit;
  
  %* Correlation table;
  proc corr data=&insdn out=tmp noprint;
    var &dep_var &varlist_redu;
  run;
  proc sql; 
    create table corrmatrix (drop=_type_) as 
    select * from tmp where _type_ = 'CORR' and _name_ ^= "&dep_var" order by abs(&dep_var) desc;
  quit;

  %* Cycle through to eliminate variable wiht correlations that are too high;
  %let breakdelcorr=0;
  %do %while(&breakdelcorr=0);
    data corrmatrix (drop=i); set corrmatrix;
      varnum = _n_;
      cnt = 0;
      array vars{*} &prefix.: ;
      do i = 1 to dim(vars);
        if &maxcorr < abs(vars{i}) then cnt = cnt + 1;
      end;
    run;

	 proc sql noprint; 
		 select count(*) into: vcnt from corrmatrix where cnt>1;
	 quit;

	 %if &vcnt = 0 %then %let breakdelcorr=1;  %* No more variables to drop;
	 %else %do;
    	%let ndropvar=0;
    	data _null_; set corrmatrix (where=(cnt>1) obs=1) end=eof;
      	call symputx("keepvar",_name_);
      	if eof then call symputx("ndropvar",1);
    	run;

    	data _null_; set corrmatrix (where=(abs(&keepvar) > &maxcorr and _name_ ^= "&keepvar")) end=lastob;
      	call symputx('dropvar'||left(_n_),_name_);
      	if lastob then call symputx('ndrop',_n_);
    	run;

    	%if &ndropvar=0 %then %let breakdelcorr=1;  %* No more variables to drop;
    	%else %do;  %* More variables to drop;
      	data corrmatrix; set corrmatrix (drop= %do i = 1 %to &ndrop; &&dropvar&i %end; );
        	%do i = 1 %to &ndrop;
          	if upcase(_name_) = upcase("&&dropvar&i") then delete;
        	%end;
      	run;
	  		data out.CE3_Var_Redu; set out.CE3_Var_Redu;
        		%do i = 1 %to &ndrop;
	          	if upcase(variable) = upcase("&&dropvar&i") then do;
            		drop_corr = 'Y';
          		end;
        		%end;
       	run;
		%end;
    %end;
  %end;
  proc datasets nolist; delete tmp corrmatrix; quit; run;
%mend;

**********************************************************************************************************************;
*** 4.Model selection and tuning ***;
**********************************************************************************************************************;
** Master macro for binary dependent variable **;
%macro Model_Val_Logistic(insdnT, insdnV, regvlist, sle=0.05, sls=.05,
			metric_sdnout=out.CE4_metric_sdnout);

	%* Set up independent variable lists;
	%if &excludelist =  %then %let nexcludes = 0;
	  %else %let nexcludes = %sysfunc(countw(&excludelist));  /* get number of variables in exclude list */
	%if &nexcludes ^= 0 %then %do;
		%let regvlist=%StringMinus(&regvlist, &excludelist);  /* Remove excluded variables from variable list */
	%end;

	%if &startlist =  %then %let nstarts = 0;
	  %else %let nstarts = %sysfunc(countw(&startlist));  /* get number of variables in start list */
	%if &nstarts ^= 0 %then %do;
		%if nexcludes ^= 0 %then %do;  /* Make sure there are no excluded variables in the start list*/
			%let startlist=%StringMinus(&startlist, &excludelist);
		%end;
		%let nstarts = %sysfunc(countw(&startlist));  /* re-get number of variables in start list */
		%if &nstarts ^= 0 %then %do;
			%let regvlist=&startlist %StringMinus(&regvlist, &startlist);
		%end;
	%end;

	%if &includelist =  %then %let nincludes = 0;
	  %else %let nincludes = %sysfunc(countw(&includelist));  /* get number of variables in include list */
	%if &nincludes ^= 0 %then %do; 
		%if nexcludes ^= 0 %then %do;  /* Make sure there are no excluded variables in the include list*/
			%let includelist=%StringMinus(&includelist, &excludelist);
		%end;
		%let nincludes = %sysfunc(countw(&includelist));  /* re-get number of variables in include list */
		%if &nincludes ^= 0 %then %do;
			%let regvlist=&includelist %StringMinus(&regvlist, &includelist);
		%end;
	%end;

	ods listing close;
	
	%** Build forward model **;
	ods output  ParameterEstimates=Parm1;
	proc logistic data=&insdnT desc namelen=32;
		model &dep_var = &regvlist /selection=forward lackfit rsq details slentry=&sle
		%if &nincludes ne 0 %then %do; include=&nincludes %end;
		%if &nstarts ne 0 %then %do; start=&nstarts %end;
		;
	run;

	%** Build backward model **;
	ods output  ParameterEstimates=Parm2;
	proc logistic data=&insdnT desc namelen=32;
		model &dep_var = &regvlist /selection=backward lackfit rsq details slstay=&sls
		%if &nincludes ne 0 %then %do; include=&nincludes %end;
		%if &nstarts ne 0 %then %do; start=&nstarts %end;
		;
	run;
	ods listing;

	%** Extract variables **;
	proc sql noprint; 
		select max(step) into: max1 from Parm1;
		select max(step) into: max2 from Parm2;
	quit;
	data parm1; set parm1 (where=(step=&max1)); run;
	data parm2; set parm2 (where=(step=&max2)); run;

	proc sql noprint;
		create table varlisttmp as
		select variable
		from parm1
		union
		select variable 
		from parm2;
		select variable into: varlisttmp separated by ' ' from varlisttmp where variable ^= 'Intercept';
		select count(*) into: nvars from varlisttmp where variable ^= 'Intercept';
	quit;
	%put &varlisttmp;

	%ctl_Stats_Log(&insdnT, &varlisttmp, &nvars);
	data &metric_sdnout;
		length file $4.;
		set stats;
		file = 'Dev';
	run;
	%ctl_Stats_Log(&insdnV, &varlisttmp, &nvars);
	data &metric_sdnout;
		set &metric_sdnout stats (in=a);
		if a then file = 'Val';
	run;

	%** Clean up;
	proc datasets library=work nolist;
		delete varlisttmp globalfit full vars tmp1 tmp2 parms vif stats;
	run;
	quit;

%mend;

** Master macro for continuous dependent variable **;
%macro Model_Val_Reg(insdnT, insdnV, regvlist, sle=0.05, sls=.05,
			metric_sdnout=out.CE4_metric_sdnout);

	%* Set up independent variable lists;
	%if &excludelist =  %then %let nexcludes = 0;
	  %else %let nexcludes = %sysfunc(countw(&excludelist));  /* get number of variables in exclude list */
	%if &nexcludes ^= 0 %then %do;
		%let regvlist=%StringMinus(&regvlist, &excludelist);  /* Remove excluded variables from variable list */
	%end;

	%if &startlist =  %then %let nstarts = 0;
	  %else %let nstarts = %sysfunc(countw(&startlist));  /* get number of variables in start list */
	%if &nstarts ^= 0 %then %do;
		%if nexcludes ^= 0 %then %do;  /* Make sure there are no excluded variables in the start list*/
			%let startlist=%StringMinus(&startlist, &excludelist);
		%end;
		%let nstarts = %sysfunc(countw(&startlist));  /* re-get number of variables in start list */
		%if &nstarts ^= 0 %then %do;
			%let regvlist=&startlist %StringMinus(&regvlist, &startlist);
		%end;
	%end;

	%if &includelist =  %then %let nincludes = 0;
	  %else %let nincludes = %sysfunc(countw(&includelist));  /* get number of variables in include list */
	%if &nincludes ^= 0 %then %do; 
		%if nexcludes ^= 0 %then %do;  /* Make sure there are no excluded variables in the include list*/
			%let includelist=%StringMinus(&includelist, &excludelist);
		%end;
		%let nincludes = %sysfunc(countw(&includelist));  /* re-get number of variables in include list */
		%if &nincludes ^= 0 %then %do;
			%let regvlist=&includelist %StringMinus(&regvlist, &includelist);
		%end;
	%end;

	ods listing close;
	ods graphics on;
	%** Build forward model **;
	ods output SelParmEst=Parm1;
	proc reg data=&insdnT plots(maxpoints=none);
		model &dep_var = &regvlist /selection=forward details slentry=&sle
		%if &nincludes ne 0 %then %do; include=&nincludes %end;
		%if &nstarts ne 0 %then %do; start=&nstarts %end;
		;
	run;
	ods listing;
	ods graphics off;
	quit;

	ods listing close;
	ods graphics on;
	%** Build backward model **;
	ods output  SelParmEst=Parm2;
	proc reg data=&insdnT plots(maxpoints=none);
		model &dep_var = &regvlist /selection=backward details slstay=&sls
		%if &nincludes ne 0 %then %do; include=&nincludes %end;
		%if &nstarts ne 0 %then %do; start=&nstarts %end;
		;
	run;
	ods listing;
	ods graphics off;
	quit;

	%** Extract variables **;
	proc sql noprint; 
		select max(step) into: max1 from Parm1;
		select max(step) into: max2 from Parm2;
	quit;
	data parm1; set parm1 (where=(step=&max1)); run;
	data parm2; set parm2 (where=(step=&max2)); run;

	proc sql noprint;
		create table varlisttmp as
		select variable
		from parm1
		union
		select variable 
		from parm2;
		select variable into: varlisttmp separated by ' ' from varlisttmp where variable ^= 'Intercept';
		select count(*) into: nvars from varlisttmp where variable ^= 'Intercept';
	quit;
	%put &varlisttmp;

	%** Get dev model statistics **;
	%ctl_Stats_Reg(&insdnT, &varlisttmp, &nvars );
	data &metric_sdnout;
		length file $4.;
		set stats;
		file = 'Dev';
	run;
	%** Get val model statistics **;
	%ctl_Stats_Reg(&insdnV, &varlisttmp, &nvars);
	data &metric_sdnout;
		set &metric_sdnout stats (in=a);
		if a then file = 'Val';
	run;

	%* Clean up;
	proc datasets library=work nolist;
		delete varlisttmp globalfit full vars tmp1 tmp2 vif stats;
	run;
	quit;

%mend;

**  StringMinus substracts one string from another  **;
%macro StringMinus(string1, string2);
	%local count word StringOut;
	%let count=1;
	%let word=%scan(&string1,&count,%str( ));
	%do %while(&word ne);
		%let count=%eval(&count+1);
		%if %sysfunc(indexw(%upcase(&string2), %upcase(&word)))=0 %then %do;
			%let stringout=&stringout &word;
		%end;
		%let word=%scan(&string1,&count,%str(  ));
	%end;
	&StringOut
%mend;

**  Remember_nobs_ngoods_nbads gets record counts  **;
%macro Remember_nobs_ngoods_nbads(insdn,bads=1,goods=0);
	%global vnobs;
	%global vngoods;
	%global vnbads;
	%if &weight= %then %do;
		proc sql noprint; 
			select count(*) into: vnobs from &insdn;
			select count(*) into: vnbads from &insdn where &dep_var in (&bads);
			select count(*) into: vngoods from &insdn where &dep_var in (&goods);
		quit;
	%end;
	%else %do;
		proc sql noprint;
			select sum(&weight) into: vnobs from &insdn;
			select sum(&weight) into: vnbads from &insdn where &dep_var in (&bads);
			select sum(&weight) into: vngoods from &insdn where &dep_var in (&goods);
		quit;
	%end;
%mend;

**  Control program for binary stats  **;
%macro ctl_Stats_Log(insdn, varlisttmp, nvars);

	%** Get counts for validation file **;
	%Remember_nobs_ngoods_nbads(&insdn);

	%** Get full model statistics **;
	%Stats_Log(&insdn, &varlisttmp, Full Model, int=N, vif=Y);
	data full; set globalfit; run;

	%** Get model statistics without intercept **;
	%Stats_Log(&insdn, &varlisttmp, Intercept, int=Y, vif=N);
	data vars; set globalfit; run;

	%** Get model statistics without specific variables **;
	%do I = 1 %to &nvars;
		%let var = %scan(&varlisttmp,&I);
		%let tmplist=%StringMinus(&varlisttmp, &var);  /* Remove single variable */
		%Stats_Log(&insdn, &tmplist, &var, int=N, vif=N);
		data vars; set vars globalfit; run;
	%end;

	%** Create macro variables for full models metrics **;
	data _null_; set full;
		call symputx('AIC',AIC);
		call symputx('SC',SC);
		call symputx('LogL2',LogL2);
		call symputx('Rsquare',Rsquare);
		call symputx('SomersD',SomersD);
		call symputx('Gamma',Gamma);
		call symputx('TauA',TauA);
		call symputx('c',c);
		call symputx('Concord',Concord);
		call symputx('Discon',Discon);
		call symputx('LackFit',LackFit);
		call symputx('Lift',Lift_index);
		call symputx('infv',infv);
		call symputx('ks',ks);
	run;
	proc sql noprint;
		select sum(abs(standardizedest)) into: cum from parms where variable ^= 'Intercept';
	quit;

	%** Subtract full model metrics from sub-model metrics;
	%** This shows the impact of each variable;
	data tmp1; set vars;
		AIC = AIC - &AIC;
		SC = SC - &SC;
		LogL2 = LogL2 - &LogL2;
		Rsquare = Rsquare - &Rsquare;
		SomersD = SomersD - &SomersD;
		Gamma = Gamma - &Gamma;
		TauA = TauA - &TauA;
		c = c - &c;
		Concord = Concord - &Concord;
		Discon = Discon - &Discon;
		LackFit = LackFit - &LackFit;
		lift_Index = Lift_Index - &Lift;
		infv = infv - &infv;
		ks = ks - &ks;
	run;

	%** Combine datasets **;
	data tmp2; set full (in=a) tmp1;
		if a then sord = 1;
			else if variable = 'Intercept' then sord = 2;
			else sord = 3;	
	run;
	proc sql;
		create table stats as
		select a.Variable,
				 b.Label,
				 b.Estimate,
				 b.StandardizedEst as StdEst,
				 b.StdErr,
				 b.WaldChiSq,
				 b.ProbChiSq,
				 c.VarianceInflation as VIF,
				 abs(b.StandardizedEst)/&cum as RelImp,
				 a.AIC,
				 a.SC,
				 a.LogL2,
				 a.Rsquare,
				 a.SomersD,
				 a.Gamma,
				 a.TauA,
				 a.c,
				 a.Concord,
				 a.Discon,
				 a.LackFit,
				 a.Lift_Index,
				 a.infv,
				 a.ks,
				 a.sord
		from tmp2 a
		left join parms b on a.variable = b.variable
		left join vif c on a.variable = c.variable
		order by a.sord, RelImp desc;
	quit;

	%** Clean up;
	proc datasets library=work nolist;
		delete globalfit full vars tmp1 tmp2 parms vif;
	run;
	quit;

%mend;

**  Stats_Log gets statistics for binary models  **;
%macro Stats_Log(insdn, varlist, var, int=N, vif=N ,bads=1,goods=0);

	%* Extract standard model statistics;
	ods listing close;
	ods output  FitStatistics= FitStatistics;
	ods output  RSquare= RSquare;
	ods output  Association= Association;
	ods output  LackFitChiSq= LackFitChiSq;
	%if &vif=Y %then %do;
		ods output  ParameterEstimates=Parms;
	%end;
	proc logistic data=&insdn desc namelen=32;
		%if &weight ne %then %do; weight &weight; %end;
		model &dep_var = &varlist / stb rsq lackfit parmlabel
		%if &int=Y %then %do; noint %end;
		;
		output out=scoretmp p=pred;
	run;
	ods listing;

	data globalfit (drop=Criterion label1 label2);
		length variable $32.;
		variable = "&var";
		merge
			%if &int=N %then %do;
				FitStatistics (where=(Criterion='AIC') rename=(InterceptAndCovariates=AIC) keep=Criterion InterceptAndCovariates)
				FitStatistics (where=(Criterion='SC') rename=(InterceptAndCovariates=SC) keep=Criterion InterceptAndCovariates) 
				FitStatistics (where=(Criterion='-2 Log L') rename=(InterceptAndCovariates=LogL2) keep=Criterion InterceptAndCovariates)
			%end;
			%else %do;
				FitStatistics (where=(Criterion='AIC') rename=(WithCovariates=AIC) keep=Criterion WithCovariates)
				FitStatistics (where=(Criterion='SC') rename=(WithCovariates=SC) keep=Criterion WithCovariates) 
				FitStatistics (where=(Criterion='-2 Log L') rename=(WithCovariates=LogL2) keep=Criterion WithCovariates)
			%end;
				RSquare (rename=(nValue1=Rsquare) keep=nValue1)
				Association (where=(label2="Somers' D") rename=(nvalue2=SomersD) keep=label2 nvalue2)
				Association (where=(label2='Gamma') rename=(nvalue2=Gamma) keep=label2 nvalue2)
				Association (where=(label2='Tau-a') rename=(nvalue2=TauA) keep=label2 nvalue2)
				Association (where=(label2='c') rename=(nvalue2=c) keep=label2 nvalue2)
				Association (where=(label1='Percent Concordant') rename=(nvalue1=Concord) keep=label1 nvalue1)
				Association (where=(label1='Percent Discordant') rename=(nvalue1=Discon) keep=label1 nvalue1)
				LackFitChiSq (rename=(ProbChiSq=LackFit) keep=ProbChiSq);
		label AIC=" ";
		label SC=" ";
		label LogL2=" ";
		label LackFit=" ";
	run;

	%* Assign groups;
	data tmp0; set scoretmp (keep=&dep_var pred &weight);
		length group $4.;
		bad = (&dep_var = &bads);
		good = (&dep_var = &goods);
		if &dep_var = &bads then group = 'bad';
			else if &dep_var = &goods then group = 'good';
		%if &weight ^=  %then %do; xxx=round(&weight); %end;
	run;

	%* Create deciles based on predicted value;
	%if &weight = %then %do;
		proc rank data=tmp0 out=tmp1 groups=10;
			var pred;
			ranks rank;
		run;
	%end;
	%else %do;
		proc sort data=tmp0 out=tmp1; by pred; run;
		data tmp1; set tmp1;
			retain cum;
			cum + &weight;
			rank = (floor(cum*10/(&vnobs+1)));
		run;
	%end;

	%* Summarize file by rank;
	proc summary data=tmp1;
		var &dep_var pred good bad;
		class rank;
		%if &weight ne %then %do; weight &weight; %end;
		%if &weight ne %then %do; output out=tmp2 (drop= _freq_) sumwgt=count %end;
		%else %do; output out=tmp2 (rename= _freq_=count) %end;
			mean(&dep_var)=mean_dv mean(pred)=mean_pred sum(bad)=badcnt sum(good)=goodcnt;
	run;
	data _null_; set tmp2 (where=(_type_=0));
		call symputx('Norm',mean_dv);
	run;

	%* Perform Kolmogorov-Smirnov test;
	proc npar1way data=tmp1 edf noprint;
		class &dep_var;
		var pred;
		%if &weight ne %then %do; freq xxx; %end;
		output edf out=kolsmir;
	run;
	proc sql noprint; select round(_d_,.000001) into: kstmp from kolsmir; quit;

	%* Create custom statistics;
	data tmp3 (keep=Lift_Index infv ks);
		set tmp2 (where=(_type_=1)) end=eof;
		Lift_Index = mean_dv / &norm * 100;
		retain infv cumbad cumgood;
		badratio=badcnt*(1.0/&vnbads);
		goodratio=goodcnt*(1.0/&vngoods);
		if _n_=1 then do;
			infv=0.0;
			cumgood=0.0;
			cumbad=0.0;
		end;

		%* Information_value;
		if badratio+goodratio>=0.000000001 then do;
			if min(badratio,goodratio)<=0.000000001 then do;
				infv=infv+count*max(badratio,goodratio)/&vnobs;
			end;
			else do;
				infv=infv+(badratio-goodratio)*log(badratio/goodratio);
			end;
		end;
		cumbad = cumbad + badratio;
		cumgood = cumgood + goodratio;

		%* K-S value;
		ks=&kstmp;

		if eof;
	run;

	data globalfit; merge globalfit tmp3; run;

	%* Variance inflation factor;
	%if &vif=Y %then %do;
		data scoretmp; set scoretmp;
			weight=pred*(1-pred);
		run;

		ods listing close;
		ods graphics on;
		ods output ParameterEstimates=VIF;
		proc reg data=scoretmp plots=none;
			weight weight;
			model &dep_var = &varlist / vif;
		run;
		ods graphics off;
		ods listing;
	%end;

	%* Clean up;
	proc datasets library=work nolist;
		delete FitStatistics RSquare Association LackFitChiSq
				 scoretmp tmp0 tmp1 tmp2 tmp3 kolsmir;
	run;
	quit;

%mend;

**  Control program for continuous stats  **;
%macro ctl_Stats_Reg(insdn, varlisttmp, nvars);

	%global vnobs;
	%if &weight= %then %do;
		proc sql noprint; 
			select count(*) into: vnobs from &insdn;
		quit;
	%end;
	%else %do;
		proc sql noprint;
			select sum(&weight) into: vnobs from &insdn;
		quit;
	%end;

	%** Get full model statistics **;
	%Stats_Reg(&insdnV, &varlisttmp, Full Model, int=N, vif=Y);
	data full; set globalfit; run;

	%** Get model statistics without intercept **;
	%Stats_Reg(&insdnV, &varlisttmp, Intercept, int=Y, vif=N);
	data vars; set globalfit; run;

	%** Get model statistics without specific variables **;
	%do I = 1 %to &nvars;
		%let var = %scan(&varlisttmp,&I);
		%let tmplist=%StringMinus(&varlisttmp, &var);  /* Remove single variable */
		%Stats_Reg(&insdnV, &tmplist, &var, int=N, vif=N);
		data vars; set vars globalfit; run;
	%end;


	%** Create macro variables for full models metrics **;
	data _null_; set full;
		call symputx('Rsquare',Rsquare);
		call symputx('AdjRsq',AdjRsq);
		call symputx('RMSE',RMSE);
		call symputx('CoeffVar',CoeffVar);
		call symputx('AIC',AIC);
		call symputx('SBC',SBC);
		call symputx('PC',PC);
		call symputx('JP',JP);
		call symputx('Lift_Index',Lift_Index);
		call symputx('gini',gini);
		call symputx('infv',infv);
		call symputx('ks',ks);
	run;
	proc sql noprint;
		select sum(abs(standardizedest)) into: cum from vif where variable ^= 'Intercept';
	quit;

	%** Subtract full model metrics from sub-model metrics;
	%** This shows the impact of each variable;
	data tmp1; set vars;
		Rsquare = Rsquare - &Rsquare;
		AdjRsq = AdjRsq - &AdjRsq;
		RMSE = RMSE - &RMSE;
		CoeffVar = CoeffVar - &CoeffVar;
		AIC = AIC - &AIC;
		SBC = SBC - &SBC;
		PC = PC - &PC;
		JP = JP - &JP;
		Lift_Index = Lift_Index - &Lift_Index;
		gini = gini - &gini;
		infv = infv - &infv;
		ks = ks - &ks;
	run;

	%** Combine datasets **;
	data tmp2; set full (in=a) tmp1;
		if a then sord = 1;
			else if variable = 'Intercept' then sord = 2;
			else sord = 3;	
	run;
	proc sql;
		create table stats as
		select a.Variable,
				 b.Label,
				 b.Estimate,
				 b.StandardizedEst as StdEst,
				 b.StdErr,
				 b.tValue,
				 b.Probt,
				 b.VarianceInflation as VIF,
				 abs(b.StandardizedEst)/&cum as RelImp,
				 a.Rsquare,
				 a.AdjRsq,
				 a.RMSE,
				 a.CoeffVar,
				 a.AIC,
				 a.SBC,
				 a.PC,
				 a.JP,
				 a.Lift_Index,
				 a.gini,
				 a.infv,
				 a.ks,
				 a.sord
		from tmp2 a
		left join vif b on a.variable = b.variable
		order by a.sord, RelImp desc;
	quit;

	%** Clean up;
	proc datasets library=work nolist;
		delete globalfit full vars tmp1 tmp2 vif;
	run;
	quit;

%mend;

**  Stats_Reg gets statistics for continuous models  **;
%macro Stats_Reg(insdn, varlist, var, int=N, vif=N);

	%* Extract standard model statistics;
	ods listing close;
	ods graphics on;
	ods output  FitStatistics=FitStatistics;
	%if &vif=Y %then %do;
		options label;
		ods output ParameterEstimates=VIF;
	%end;
	proc reg data=&insdn outest=outest plots=none;
		%if &weight ne %then %do; weight &weight; %end;
		model &dep_var = &varlist / stb aic sbc jp pc
		%if &vif=Y %then %do; vif %end;
		%if &int=Y %then %do; noint %end;
		;
		output out=scoretmp p=pred;
	run;
	ods graphics off;
	ods listing;
	quit;

	data globalfit (drop=label1 label2);
		length variable $32.;
		variable = "&var";
		merge FitStatistics (where=(label2="R-Square") rename=(nvalue2=Rsquare) keep=label2 nvalue2)
				FitStatistics (where=(label2="Adj R-Sq") rename=(nvalue2=AdjRsq) keep=label2 nvalue2)
				FitStatistics (where=(label1="Root MSE") rename=(nvalue1=RMSE) keep=label1 nvalue1)
				FitStatistics (where=(label1="Coeff Var") rename=(nvalue1=CoeffVar) keep=label1 nvalue1)
				outest (keep=_AIC_ _SBC_ _PC_ _JP_ rename=(_AIC_=AIC _SBC_=SBC _PC_=PC _JP_=JP));
	run;

	%* Create deciles based on predicted value;
	%if &weight = %then %do;
		proc rank data=scoretmp (keep=&dep_var pred &weight) out=tmp1 groups=10;
			var pred;
			ranks rank;
		run;
	%end;
	%else %do;
		proc sort data=scoretmp (keep=&dep_var pred &weight) out=tmp1; by pred; run;
		data tmp1; set tmp1;
			retain cum;
			cum + &weight;
			rank = (floor(cum*10/(&vnobs+1)));
		run;
	%end;

	%* Summarize file by rank;
	proc summary data=tmp1;
		var &dep_var pred;
		class rank;
		%if &weight ne %then %do; weight &weight; %end;
		%if &weight ne %then %do; output out=tmp2 (drop= _freq_) sumwgt=count %end;
		%else %do; output out=tmp2 (rename= _freq_=count) %end;
			mean(&dep_var)=mean_dv mean(pred)=mean_pred;
	run;

	%* Create custom statistics;
	proc sql noprint;
		select mean_dv into: norm   from tmp2 (where=(_type_=0));
		select count into: ntotal from tmp2 (where=(_type_=0));
		select sum(count*mean_dv) into: totalresp from tmp2 (where=(_type_=1));
		select sum(count*(count*mean_dv/&totalresp)) into: totalbad from tmp2 (where=(_type_=1));
	quit;

	data tmp3 (keep=lift_index gini infv ks);
		set tmp2 (where=(_type_=1)) end=eof;
		Lift_Index = mean_dv / &norm * 100;
		retain cumresp cumtotal;
		if _n_ = 1 then do;
			cumresp = 0;
			cumtotal = 0;
		end;
		cumresp = cumresp + mean_dv * count;
		cumpct_resp = cumresp / &totalresp;
		cumtotal = cumtotal + count;
		cumpct_freq = cumtotal / &ntotal;
		cumavg = cumresp / cumtotal;
		cumindex = cumavg / &norm * 100;
		badrate = count * mean_dv / &totalresp;
		goodrate = (1 - badrate) * count / (&ntotal - &totalbad);

		retain gini infv ks cumbad cumgood 0;
		prev_dv = lag(cumpct_resp);
		prev_pct = lag(cumpct_freq);
		if _n_ ne 1 then gini = gini+2*((cumpct_freq+prev_pct)/2-(cumpct_resp+prev_dv)/2)*(cumpct_freq-prev_pct);
		infv = infv+(badrate-goodrate)*log(badrate/goodrate);
		cumbad = cumbad+badrate;
		cumgood = cumgood+goodrate;
		ks = max(ks,abs(cumbad-cumgood));
		if eof;
	run;

	data globalfit; merge globalfit tmp3; run;

	%* Clean up;
	proc datasets library=work nolist;
		delete FitStatistics outest scoretmp tmp1 tmp2 tmp3;
	run;
	quit;

%mend;

** Variable Tuning and Selection Based on model metric output **;
%macro Vars_tune(dt1=, dt2=, out_txt=&path_output.CE4_Varlist_Final.txt);

	%* Check criteria and reset if necessary;
	%if %upcase(&binary_dv) = Y %then %do;
		%let c=AIC SC LOGL2 RSQUARE SOMERSD GAMMA TAUA C CONCORD DISCON LACKFIT LIFT_INDEX INFV KS;
		%if %sysfunc(indexw(&c,%upcase(&criteria)))=0 %then %do;
			%let criteria = c;
			%put Criteria changed to c;
		%end;
	%end;
	%else %do;
		%let c=AIC SBC JP RSQUARE ADJRSQ RMSE COEFFVAR PC LIFT_INDEX GINI INFV KS;
		%if %sysfunc(indexw(&c,%upcase(&criteria)))=0 %then %do;
			%let criteria = AdjRsq;
			%put Criteria changed to AdjRsq;
		%end;
	%end;

	%let neg=RSQUARE SOMERSD GAMMA TAUA C CONCORD INFV TVR KS DIVERG ADJRSQ LIFT_INDEX GINI;
	%if &threshold ^= 0 and %index(&neg,%upcase(&criteria))>0 %then %do;
		%let threshold = %sysevalf(-1*&threshold);
	%end;

	proc sql;
		create table tempfinal as
		select variable 
		from &dt1
		where sord = 3 and &criteria <= &threshold
		  and RelImp >= &MinImp
		  and file = 'Val'
		&SQL_join
		select variable 
		from &dt2
		where sord = 3 and &criteria <= &threshold
		  and RelImp >= &MinImp
		  and file = 'Val';
	quit;

	data _NULL_;
		set tempfinal end=eof;
		FILE  "&out_txt"  LRECL=256;
		if _N_=1 then do;
			PUT ' ';
			PUT "%"@;
			PUT "LET Varlist_Final=";
		end;
		PUT variable @;
		if eof then PUT ';';
	run;

	proc sort data=tempfinal; by variable;
	data out.ce4_variables;
		merge out.ce4_variables tempfinal (in=a);
		by variable;
		if a then final = 'Y';
	run;

	proc datasets library=work nolist;
		delete tempfinal;
	quit;
	run;

%mend;

** Create Graph report **;
%macro Grafing (inds,varlist,cats);

	%** Get variables and unique counts **;
	proc sql;
		create table cnts as
		select 		   
			%let i=1;
			%let v=%scan(&varlist, &i,' ');
			%do %while (&v^=);
				count(distinct(&v)) as &v,
				%let i=%eval(&i+1);
				%let v=%scan(&varlist, &i,' ');
			%end;
			"dummy" as dummy
		from &inds;
	quit;
	proc transpose data=cnts (drop=dummy) out=cnts (rename=col1=uniq) name=variable; run;

	%** Create macro variables to control processing **;
	data _null_; set cnts end=eof;
		if eof then call symputx("VARCNT",_N_);
		call symputx("IV"|| trim(left(put(_N_,4.)))  ,variable);
		call symputx("uniq"|| trim(left(put(_N_,4.)))  ,uniq);
	run;

	goptions reset=all device=pdf display gunit=pct border ftext= htitle=8 htext=3;
	ods listing close;
	ods pdf file="&path_output.CE4_Graphs.pdf" style=sasweb startpage=no;
	%** Loop through for each variable **;
	%do _I_ = 1 %to &varcnt;
		%if &_I_ ^= 1 %then %do; ods pdf startpage=now; %end;
		%if &&uniq&_I_ <= &cats %then %do;
			proc summary data=&inds nway missing;
				var &dep_var;
				class &&iv&_I_;
				output out=graf (drop=_type_ rename=_freq_=count) mean=mean;
			run;

			proc print data=graf noobs;
				title "Variable: &&iv&_I_";
				var &&iv&_I_ count mean;
				format count comma8.;
				format mean 8.4;
			run;
			
			proc sgplot data=graf;
				vbar &&iv&_I_ / response=mean nostatlabel /*nooutline*/ fillattrs=(color="lightblue") /*transparency=.5*/;
				vline &&iv&_I_ / response=count nostatlabel y2axis lineattrs=(color="darkblue" thickness=2);
				label count="# of Customers";
				label mean="Mean &dep_var";
				label &&iv&_I_="&&iv&_I_";
				keylegend / location = outside
				position = top
				noborder
				title = "&&iv&_I_";
				format count comma8.;
			run;
		%end;
		%else %do;
			proc means data=&inds p1 p99 NOPRINT;
     			var &&&iv&_I_;
     			output out=tmp p1=var_p1 p99=var_p99  / noinherit;
  			run;
			data _NULL_;
				set tmp;
				range = (var_p99 - var_p1)/&cats;
				call symputx('cut_lo',var_p1);
				call symputx('cut_hi',var_p99);
				call symputx('range',range);
			run;

			data tmp;
				set &inds (keep=&dep_var &&iv&_I_);
				if &&iv&_I_ < &cut_lo then bin = 1;
				else if &&iv&_I_ >= &cut_hi then bin = input("&cats",best10.);
				else do;
					do k = 1 to &cats;
						if &&iv&_I_ >= &cut_lo+(k-1)*&range and &&iv&_I_< &cut_lo+k*&range then	bin=k;
					end;
				end;
			run;

			proc summary data=tmp nway missing;
				var &&iv&_I_ &dep_var;
				class bin;
				output out=graf (drop=_type_ rename=_freq_=count) min(&&iv&_I_)=lo max(&&iv&_I_)=hi mean(&dep_var)=mean;
			run;

			proc print data=graf noobs;
				title "Variable: &&iv&_I_";
				var bin count lo hi mean;
				format count comma8.;
				format lo hi mean 8.4;
			run;

			proc sgplot data=graf;
				vbar lo / response=mean nostatlabel fillattrs=(color="lightblue");
				vline lo / response=count nostatlabel y2axis lineattrs=(color="darkblue" thickness=2);
				label count="# of Customers";
				label mean="Mean &dep_var";
				label lo="&&iv&_I_";
				keylegend / location = outside
				position = top
				noborder
				title = "&&iv&_I_";
				format count comma8.;
			run;

			proc datasets nolist;
				delete tmp;
			run;
		%end;
		
	%end;
	ods pdf close;
	ods listing;

	proc datasets nolist;
		delete cnts graf;
	run;

%mend;

** Master macro for Model selection and tuning **;
%macro CE_Model_Val(insdn, varlist);

	data mod val;
		set &insdn (where=(mod_val_test ^= 3));
		keep &Dep_var mod_val_test &varlist &weight;
		if mod_val_test=1 then output mod;
			else if mod_val_test=2 then output val;
	run;

	%* Start variable table;
	data out.ce4_variables;
		length variable $32;
      %Let I = 1;
      %Let var = %scan(&varlist,&I);
      %Do %while(&var ne );
            variable="&var";
            output;
         %Let I = %eval(&I + 1);
         %Let var = %scan(&varlist,&I);
      %end;
	run;
	proc sort data=out.ce4_variables; by variable; run;

	%if %upcase(&Binary_dv) = Y %then %do;
		%Model_Val_Logistic(mod, val, &varlist,
				sle=&sel_alpha,
				sls=&sel_alpha,
				metric_sdnout=out.CE4_Model_Metric_Mod);

		proc sort data=parm1 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		proc sort data=parm2 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		data out.ce4_variables;
			merge out.ce4_variables parm1 (in=a) parm2 (in=b);
			by variable;
			if a then Tforward = 'Y';
			if b then Tbackward = 'Y';
		run;
		%Model_Val_Logistic(val, mod, &varlist,
				sle=&sel_alpha,
				sls=&sel_alpha,
				metric_sdnout=out.CE4_Model_Metric_Val);

		proc sort data=parm1 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		proc sort data=parm2 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		data out.ce4_variables;
			merge out.ce4_variables parm1 (in=a) parm2 (in=b);
			by variable;
			if a then Vforward = 'Y';
			if b then Vbackward = 'Y';
		run;

		%*model selection and tuning;
		%Vars_tune(dt1=out.CE4_Model_Metric_Mod,
				dt2=out.CE4_Model_Metric_Val,
				out_txt=&path_output.CE4_Varlist_Final.txt);

		%* Report;
		ods listing close;
		ods Tagsets.ExcelxP body="&path_output.CE4_Model_report.xls" style=sasweb;

		ods tagsets.excelxp options(sheet_name="Variables");
		proc print data=out.ce4_variables noobs; run;

		ods tagsets.excelxp options(sheet_name="Model" embedded_titles="yes");
		proc print data=out.CE4_Model_Metric_Mod noobs;
			title "Diagnostic statistics: Model Portion";
			var file Variable Label;
			var Estimate  / style={tagattr='format:0.000000'};
			var StdEst StdErr WaldChiSq ProbChiSq VIF RelImp  / style={tagattr='format:0.0000'};
			var AIC SC LogL2  / style={tagattr='format:0.00'};
			var Rsquare SomersD Gamma TauA c Concord Discon LackFit / style={tagattr='format:0.0000'};
			var Lift_Index / style={tagattr='format:0.00'};
			var infv ks  / style={tagattr='format:0.0000'};
		run;
		title ;

		ods tagsets.excelxp options(sheet_name="Validation" embedded_titles="yes");
		proc print data=out.CE4_Model_Metric_Val noobs;
			title "Diagnostic statistics: Validation Portion";
			var file Variable Label;
			var Estimate  / style={tagattr='format:0.000000'};
			var StdEst StdErr WaldChiSq ProbChiSq VIF RelImp  / style={tagattr='format:0.0000'};
			var AIC SC LogL2  / style={tagattr='format:0.00'};
			var Rsquare SomersD Gamma TauA c Concord Discon LackFit / style={tagattr='format:0.0000'};
			var Lift_Index / style={tagattr='format:0.00'};
			var infv ks  / style={tagattr='format:0.0000'};
		run;
		title ;

		ods Tagsets.ExcelxP close;
		ods listing;

	%end;
	%else %if %upcase(&Binary_dv) ^= Y %then %do;
		%Model_Val_Reg(mod, val, &varlist,
				sle=&sel_alpha,
				sls=&sel_alpha,
				metric_sdnout=out.CE4_Model_Metric_Mod);

		proc sort data=parm1 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		proc sort data=parm2 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		data out.ce4_variables;
			merge out.ce4_variables parm1 (in=a) parm2 (in=b);
			by variable;
			if a then Tforward = 'Y';
			if b then Tbackward = 'Y';
		run;

		%Model_Val_Reg(val, mod, &varlist,
				sle=&sel_alpha,
				sls=&sel_alpha,
				metric_sdnout=out.CE4_Model_Metric_Val);

		proc sort data=parm1 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		proc sort data=parm2 (keep=variable where=(variable ^= 'Intercept')); by variable; run;
		data out.ce4_variables;
			merge out.ce4_variables parm1 (in=a) parm2 (in=b);
			by variable;
			if a then Vforward = 'Y';
			if b then Vbackward = 'Y';
		run;

		%*model selection and tuning;
		%Vars_tune(dt1=out.ce4_model_metric_mod,
				dt2=out.ce4_model_metric_val,
				out_txt=&path_output.CE4_Varlist_Final.txt);

		%* Report;
		ods listing close;
		ods Tagsets.ExcelxP body="&path_output.CE4_Model_report.xls" style=sasweb;

		ods tagsets.excelxp options(sheet_name="Variables");
		proc print data=out.ce4_variables noobs; run;

		ods tagsets.excelxp options(sheet_name="Model" embedded_titles="yes");
		proc print data=out.CE4_Model_Metric_Mod noobs;
			title "Diagnostic statistics: Model Portion";
			var file Variable Label;
			var Estimate  / style={tagattr='format:0.000000'};
			var StdEst StdErr tValue Probt VIF RelImp  / style={tagattr='format:0.0000'};
			var AIC  SBC JP  / style={tagattr='format:0.00'};
			var Rsquare AdjRsq RMSE CoeffVar PC / style={tagattr='format:0.0000'};
			var Lift_Index / style={tagattr='format:0.00'};
			var gini infv ks  / style={tagattr='format:0.0000'};
		run;
		title ;

		ods tagsets.excelxp options(sheet_name="Validation" embedded_titles="yes");
		proc print data=out.CE4_Model_Metric_Val noobs;
			title "Diagnostic statistics: Validation Portion";
			var file Variable Label;
			var Estimate  / style={tagattr='format:0.000000'};
			var StdEst StdErr tValue Probt VIF RelImp  / style={tagattr='format:0.0000'};
			var AIC SBC JP  / style={tagattr='format:0.00'};
			var Rsquare AdjRsq RMSE CoeffVar PC / style={tagattr='format:0.0000'};
			var Lift_Index / style={tagattr='format:0.00'};
			var gini infv ks  / style={tagattr='format:0.0000'};
		run;
		title ;

		ods Tagsets.ExcelxP close;
		ods listing;

	run;
	title ;
	%end;

	%* Do charting?;
	%if %upcase(&graph_plot)=Y %then %do;
		%inc "&path_output.CE4_Varlist_Final.txt";
		data gdat; set &insdn (where=(mod_val_test ^= 3)); run;
		%Grafing(gdat,&Varlist_Final,20);
		proc datasets library=work nolist;
			delete gdat;
		run;
		quit;
	%end;

	proc datasets library=work nolist;
		delete mod val parm1 parm2;
	run;
	quit;

%mend CE_Model_Val;

**********************************************************************************************************************;
*** 5.Final Model build and validation on test sample ***;
**********************************************************************************************************************;
** Master macro for binary dependent variable **;
%macro Model_Fin_Logistic(insdnT, insdnV, regvlist, sle=0.05, sls=.05,
			metric_sdnout=out.CE5_Model_Metric_All,
			out_txt=&path_output.CE5_Varlist_Model.txt, outds=out.CE5_scored);

	%** Build model **;
	ods listing close;
	ods output  ParameterEstimates=Parm1;
	ods output  ModelBuildingSummary= out.ce5_SelectionSummary;
	proc logistic data=&insdnT outmodel=est desc namelen=32;
		%if &weight ne %then %do; weight &weight; %end;
		model &dep_var = &regvlist /selection=&method stb lackfit rsq details
		%if &method = stepwise or &method = forward %then %do; slentry=&sle %end; 
		%if &method = stepwise or &method = backward %then %do; slstay=&sls %end; 
		;
		output out=out1 P=pscore;
	run;
	ods listing;

	%** Extract variables **;
	proc sql noprint; 
		select max(step) into: max1 from Parm1;
	quit;
	data parm1; set parm1 (where=(step=&max1)); run;

	proc sql noprint;
		select variable into: varlisttmp separated by ' ' from parm1 where variable ^= 'Intercept';
		select count(*) into: nvars from parm1 where variable ^= 'Intercept';
	quit;
	%put &varlisttmp;

	%* Score validation dataset;
	proc logistic inmodel=est;
   	score data=&insdnV out=out2;
	run;
	data &outds; set out1 (drop=_level_) out2 (drop=f_&dep_var I_&dep_var p_0 rename=(P_1=pscore)); run;

	%ctl_Stats_Log(&insdnT, &varlisttmp, &nvars);
	data &metric_sdnout;
		length file $4.;
		set stats;
		file = 'Dev';
	run;
	%ctl_Stats_Log(&insdnV, &varlisttmp, &nvars);
	data &metric_sdnout;
		set &metric_sdnout stats (in=a);
		if a then file = 'Val';
	run;

	%* Write out final variable list;
	data _NULL_;
		set parm1 (where=(variable ^= 'Intercept')) end=eof;
		FILE  "&out_txt"  LRECL=256;
		if _N_=1 then do;
			PUT ' ';
			PUT "%"@;
			PUT "LET Varlist_Final=";
		end;
		PUT variable @;
		if eof then PUT ';';
	run;

	%** Clean up;
	data out.ce5_parameterestimates; set parm1 (drop=step); run;

	proc datasets library=work nolist;
		delete parm1 stats est out1 out2;
	run;
	quit;

%mend;

** Master macro for continuous dependent variable **;
%macro Model_Fin_Reg(insdnT, insdnV, regvlist, sle=0.05, sls=.05,
			metric_sdnout=out.CE5_Model_Metric_All,
			out_txt=&path_output.CE5_Varlist_Model.txt, outds=out.CE5_scored);

	%** Build model **;
	ods listing close;
	ods graphics on;
	ods output SelParmEst=Parm1;
	ods output SelectionSummary=out.ce5_SelectionSummary;
	proc reg data=&insdnT outest=est plots(maxpoints=none);
		%if &weight ne %then %do; weight &weight; %end;
		model &dep_var = &regvlist /selection=&method stb details
		%if &method = stepwise or &method = forward %then %do; slentry=&sle %end; 
		%if &method = stepwise or &method = backward %then %do; slstay=&sls %end; 
		;
		output out=out1 P=pscore;
	run;
	ods listing;
	ods graphics off;
	quit;

	%** Extract variables **;
	proc sql noprint; 
		select max(step) into: max1 from Parm1;
	quit;
	data parm1; set parm1 (where=(step=&max1)); run;

	proc sql noprint;
		select variable into: varlisttmp separated by ' ' from parm1 where variable ^= 'Intercept';
		select count(*) into: nvars from parm1 where variable ^= 'Intercept';
	quit;
	%put &varlisttmp;

	%* Score validation dataset;
	proc score data=&insdnV score=est out=out2 type=parms;
		var &varlisttmp;
	run;
	data &outds; set out1 out2 (rename=(model1=pscore)); run;

	%ctl_Stats_Reg(&insdnT, &varlisttmp, &nvars);
	data &metric_sdnout;
		length file $4.;
		set stats;
		file = 'Dev';
	run;
	%ctl_Stats_Reg(&insdnV, &varlisttmp, &nvars);
	data &metric_sdnout;
		set &metric_sdnout stats (in=a);
		if a then file = 'Val';
	run;

	%* Write out final variable list;
	data _NULL_;
		set parm1 (where=(variable ^= 'Intercept')) end=eof;
		FILE  "&out_txt"  LRECL=256;
		if _N_=1 then do;
			PUT ' ';
			PUT "%"@;
			PUT "LET Varlist_Final=";
		end;
		PUT variable @;
		if eof then PUT ';';
	run;

	%** Clean up;
	data out.ce5_parameterestimates; set parm1 (drop=step); run;
	proc datasets library=work nolist;
		delete parm1 stats est out1 out2;
	run;
	quit;

%mend;

** Build gains table to check performance **;
%macro gains(inputfile = , score = , varlist= , title_key = );
	%* Create deciles based on predicted value;
	%if &weight = %then %do;
		proc rank data=&inputfile out=tmp1 descending groups=10;
			var &score;
			ranks decile;
		run;
	%end;
	%else %do;
		proc sql noprint;
			select sum(&weight) into: nobs from &inputfile;
		quit;
		proc sort data=&inputfile out=tmp1; by descending &score; run;
		data tmp1; set tmp1;
			retain cum;
			cum + &weight;
			decile = (floor(cum*10/(&nobs+1)));
		run;
	%end;

	%* Summarize file by decile;
	proc summary data=tmp1;
		var &dep_var &score &varlist ;
		class decile;
		%if &weight ne %then %do; weight &weight; %end;
		%if &weight ne %then %do; output out=tmp2 (drop= _freq_) sumwgt=count %end;
		%else %do; output out=tmp2 (rename= _freq_=count) %end;
			mean= min(&score)=minscore max(&score)=maxscore 
			%if &binary_dv=Y %then %do; sum(&dep_var)=numresp %end;
			;
	run;
	data out.ce5_profiles_&title_key;
		set tmp2 (drop=minscore maxscore %if &binary_dv=Y %then %do; numresp %end;);
		decile = decile + 1;
	run;

	%* Model Performance;
	%if  &binary_dv = Y %then %do;
		proc sql;
			select (numresp/count) into: avg from tmp2 where _type_=0;
			select numresp into: totresp from tmp2 where _type_=0;
		quit;
		data out.ce5_gainstable_&title_key (drop=_type_ cumresp cumcount); 
			set tmp2 (drop=&varlist where=(_type_=1));
			decile = decile + 1;
			avg_rate = numresp / count;
			lift = (avg_rate / &avg)*100;
			retain cumresp cumcount 0;
			cumresp + numresp;
			cumcount + count;
			cum_index = ((cumresp/cumcount) / &avg) * 100;
			resp_pct = numresp / &totresp;
			cum_resp_pct = cumresp / &totresp;
		run;
	%end;
	%else %do;
		proc sql;
			select &dep_var into: avg from tmp2 where _type_=0;
		quit;
		data out.ce5_gainstable_&title_key (drop=_type_); 
			set tmp2 (drop=&varlist where=(_type_=1));
			decile = decile + 1;
			lift = (&dep_var / &avg)*100;
		run;
	%end;

	%* Clean up;
	proc datasets library=work nolist;
		delete tmp1 tmp2;
	run;
	quit;
%mend;

** Profiling Control **;
%macro Fin_Profiling (inds,varlist);

	%** Initialize profiling dataset **;
	data out.CE5_profile;
		set _NULL_;
		length variable $32.;
		length label $256.;
		length category $256.;
		format count comma8.;
		format percent percent8.2;
		%if %upcase(&binary_dv) = Y %then %do;
			format Average_DV percent8.2;
		%end;
		%else %do;
			format Average_DV 12.2;
		%end;
		format index 8.0;
		length star $8.;
	run;

	%** Get overall average and count **;
	proc sql noprint;
		select count(*) into : nobs from &inds;
		select mean(&dep_var) into : overall_avg from &inds;
	quit;

	%** Get variables and unique counts **;
	proc sql;
		create table cnts as
		select 		   
			%let i=1;
			%let v=%scan(&varlist, &i,' ');
			%do %while (&v^=);
				count(distinct(&v)) as &v,
				%let i=%eval(&i+1);
				%let v=%scan(&varlist, &i,' ');
			%end;
			"dummy" as dummy
		from &inds;
	quit;
	proc transpose data=cnts (drop=dummy) out=cnts (rename=col1=uniq) name=variable; run;

	%** Create macro variables to control processing **;
	data _null_; set cnts end=eof;
		if eof then call symputx("VARCNT",_N_);
		call symputx("IV"|| trim(left(put(_N_,4.)))  ,variable);
		call symputx("uniq"|| trim(left(put(_N_,4.)))  ,uniq);
	run;
	%** Loop through for each variable **;
	%do _I_ = 1 %to &varcnt;
		%if &&uniq&_I_ <= &fin_num_category %then %do;
			%fprof1(&inds,&&iv&_I_);
		%end;
		%else %do;
			%fprof2(&inds,&&iv&_I_);
		%end;
		%fprof3(&&iv&_I_);
	%end;

	proc datasets nolist;
		delete cnts;
	run;

%mend;

** Profiling macro 1 **;
%macro fprof1(insdn,var);
	proc summary data=&insdn nway missing;
		var &dep_var;
		class &var;
		output out=prof (drop=_type_ rename=_freq_=xcount) mean=xmean;
	run;

	data prof;
		set prof;
		length xcategory $256.;
		xcategory = trim(left(put(&var,best8.)));
	run;
%mend;

** Profiling macro 2 **;
%macro fprof2(insdn, var);
	%if %upcase(&fin_equal_dist) = Y %then %do;
		proc means data=&insdn p1 p99 NOPRINT;
      	var &var;
      	output out=tmp p1=var_p1 p99=var_p99  / noinherit;
    	run;
		data _NULL_;
			set tmp;
			range = (var_p99 - var_p1)/&fin_num_category;
			call symputx('cut_lo',var_p1);
			call symputx('cut_hi',var_p99);
			call symputx('range',range);
		run;

		data tmp;
			set &insdn (keep=&dep_var &var);
			if &var < &cut_lo then bin = 1;
			else if &var >= &cut_hi then bin = input("&fin_num_category",best10.);
			else do;
				do k = 1 to &fin_num_category;
					if &var>= &cut_lo+(k-1)*&range and &var< &cut_lo+k*&range then	bin=k;
				end;
			end;
		run;
	%end;
	%else %do;
		proc rank data=&insdn (keep=&dep_var &var) out=tmp ties=High group=&fin_num_category;
			var &var;
			ranks bin;
		run;
	%end;

	proc summary data=tmp nway missing;
		var &var &dep_var;
		class bin;
		output out=prof (drop=_type_ rename=_freq_=xcount) min(&var)=lo max(&var)=hi mean(&dep_var)=xmean;
	run;

	data prof;
		set prof end=eof;
		length xcategory $256.;
		if _N_ = 1 then xcategory = "Low to " || trim(left(put(hi,best8.)));
			else if eof then xcategory = trim(left(put(lo,best8.))) || " to High";
			else xcategory = trim(left(put(lo,best8.))) || " to " || trim(left(put(hi,best8.)));
	run;

	proc datasets nolist;
		delete tmp;
	run;

%mend;

** Profiling macro 3 **;
%macro fprof3(var);

	proc sql;
		create table prof2 as
		select a.variable, c.label, b.xcategory, b.xcount, b.xmean
		from cnts a, prof b, out.ce2_vars c
		where a.variable = "&var"
		  and c.variable = "&var";
	quit;

	data prof3 (drop=xcount xmean xcategory);
		set prof2 end=eof;
		length category $256.;
		length star $8.;
		if _N_ = 1 then do;
			category = "Overall";
			Average_DV = &overall_avg;
			Count = &nobs;
			Percent = 1;
			index = 100;
			output;
		end;
		category = xcategory;
		count = xcount;
		percent = xcount / &nobs;
		Average_DV = xmean;
		index = (Average_DV / &overall_avg)*100;
		if index >= 110 then star = '* (+)';
			else if index > 100 then star = '  (+)';
			else if index <= 90 then star = '* (-)';
			else if index <= 100 then star = '  (-)';
			else star = '  (0)';
		output;
	run;

	data out.CE5_profile;
		set out.CE5_profile prof3;
	run;

	proc datasets nolist;
		delete prof prof2 prof3;
	run;

%mend;

** Master macro for Final Model build and validation **;
%macro CE_Model_Lift(insdn, varlist);

	data mod val;
		set &insdn (keep=&Dep_var mod_val_test &varlist &weight);
		if mod_val_test=3 then output val;
			else output mod;
	run;

	%if %upcase(&Binary_dv) = Y %then %do;
		%Model_Fin_Logistic(mod, val, &varlist,
				sle=&fin_alpha,
				sls=&fin_alpha,
				metric_sdnout=out.CE5_Model_Metric_All,
				out_txt=&path_output.CE5_Varlist_Model.txt,
				outds=out.CE5_scored);
	%end;
	%else %if %upcase(&Binary_dv) ^= Y %then %do;
		%Model_Fin_Reg(mod, val, &varlist,
				sle=&fin_alpha,
				sls=&fin_alpha,
				metric_sdnout=out.CE5_Model_Metric_All,
				out_txt=&path_output.CE5_Varlist_Model.txt,
				outds=out.CE5_scored);
	%end;

	%* Final variable list;
	proc sql noprint;
		select variable into: tmplst separated by ' ' from out.ce5_parameterestimates where variable ^= 'Intercept';
	quit;

	%* Performance;
	%gains(inputfile = out.CE5_scored (where=(mod_val_test ^= 3)), score = pscore,
			 varlist=&tmplst, title_key = train);

	%gains(inputfile = out.CE5_scored (where=(mod_val_test = 3)), score = pscore,
			 varlist=&tmplst, title_key = test);

	%* Profile Variables;
	%Fin_Profiling(mod,&tmplst);

	%* Report;
	ods listing close;
	ods Tagsets.ExcelxP body="&Path_output.CE5_Model_report.xls" style=sasweb;

	%* Selection Summary;
	ods tagsets.excelxp options(sheet_name="Selection");
	proc print data=out.CE5_SelectionSummary noobs; run;

	%* Parameters;
	ods tagsets.excelxp options(sheet_name="Parameters");
	proc print data=out.CE5_ParameterEstimates noobs; run;

	%* Statistics;
	%if %upcase(&Binary_dv) = Y %then %do;
		ods tagsets.excelxp options(sheet_name="Statistics");
		proc print data=out.CE5_Model_Metric_All noobs;
			var file Variable Label;
			var Estimate  / style={tagattr='format:0.000000'};
			var StdEst StdErr WaldChiSq ProbChiSq VIF RelImp  / style={tagattr='format:0.0000'};
			var AIC SC LogL2  / style={tagattr='format:0.00'};
			var Rsquare SomersD Gamma TauA c Concord Discon LackFit / style={tagattr='format:0.0000'};
			var Lift_Index / style={tagattr='format:0.00'};
			var infv ks  / style={tagattr='format:0.0000'};
		run;
	%end;
	%else %do;
		ods tagsets.excelxp options(sheet_name="Statistics");
		proc print data=out.CE5_Model_Metric_All noobs;
			var file Variable Label;
			var Estimate  / style={tagattr='format:0.000000'};
			var StdEst StdErr tValue Probt VIF RelImp  / style={tagattr='format:0.0000'};
			var AIC  SBC JP  / style={tagattr='format:0.00'};
			var Rsquare AdjRsq RMSE CoeffVar PC / style={tagattr='format:0.0000'};
			var Lift_Index / style={tagattr='format:0.00'};
			var gini infv ks  / style={tagattr='format:0.0000'};
		run;
	%end;

	%* Performance;
	ods tagsets.excelxp options(sheet_interval="table");
	ods tagsets.excelxp options(sheet_name="Performance" sheet_interval="none" embedded_titles="yes" convert_percentages="yes");
	proc print data=out.CE5_GainsTable_Train noobs;
		title 'Train';
		var decile / style={tagattr='format:0'};
		var Count  / style={tagattr='format:#,##0'};
		var &dep_var pscore / style={tagattr='format:0.0000'};
		var minscore maxscore / style={tagattr='format:0.0000'};
		%if %upcase(&Binary_dv) = Y %then %do;
			var avg_rate / style={tagattr='format:0.00%'};
			var lift cum_index / style={tagattr='format:0'};
			var resp_pct cum_resp_pct / style={tagattr='format:0.00%'};
		%end;
		%else %do;
			var lift / style={tagattr='format:0'};
		%end;
	run;
	proc print data=out.CE5_GainsTable_Test noobs;
		title 'Test';
		var decile / style={tagattr='format:0'};
		var Count  / style={tagattr='format:#,##0.'};
		var &dep_var pscore / style={tagattr='format:0.0000'};
		var minscore maxscore / style={tagattr='format:0.0000'};
		%if %upcase(&Binary_dv) = Y %then %do;
			var avg_rate / style={tagattr='format:0.00%'};
			var lift cum_index / style={tagattr='format:0'};
			var resp_pct cum_resp_pct / style={tagattr='format:0.00%'};
		%end;
		%else %do;
			var lift / style={tagattr='format:0'};
		%end;
	run;
	title ' ';

	%* Variable Validation;
	ods tagsets.excelxp options(sheet_interval="table");
	ods tagsets.excelxp options(sheet_name="Variable Validation" sheet_interval="none" embedded_titles="yes");
	proc print data=out.CE5_Profiles_Train noobs;
		title 'Train';
		var decile _type_  / style={tagattr='format:0'};
		var Count  / style={tagattr='format:#,##0'};
		var &dep_var pscore &tmplst / style={tagattr='format:0.0000'};
	run;
	proc print data=out.CE5_Profiles_Test noobs;
		title 'Test';
		var decile _type_  / style={tagattr='format:0'};
		var Count  / style={tagattr='format:#,##0'};
		var &dep_var pscore &tmplst / style={tagattr='format:0.0000'};
	run;
	title ' ';

	%* Correlations;
	proc corr data = out.ce5_scored outp=out.CE5_corr noprint; var &dep_var &tmplst; run;
	data out.CE5_corr (drop=_type_) ; set out.CE5_corr (rename=(_name_=variable) where=(_type_ = 'CORR')); run;
	ods tagsets.excelxp options(sheet_interval="table");
	ods tagsets.excelxp options(sheet_name="Correlations");
	proc print data=out.CE5_corr noobs;
		var variable;
		var &dep_var &tmplst / style={tagattr='format:0.0000'};
	run;

	%* Variable Profiles;
	proc sql noprint;
		select max(length(variable)) into: l1 from out.CE5_profile;
		select max(length(label)) into: l2 from out.CE5_profile;
		select max(length(category)) into: l3 from out.CE5_profile;
	quit;
	ods tagsets.excelxp options(sheet_name="Profiles" embedded_titles="yes" convert_percentages="yes"
		absolute_column_width="&l1,&l2,&l3,8,8,9,8,6" Frozen_Headers='Yes' Frozen_RowHeaders='1');
	proc print data=out.CE5_profile noobs;
		var variable label category;
		var count / style={tagattr='format:#,##0'};
		var percent Average_DV index star;
	run;

	ods Tagsets.ExcelxP close;
	ods listing;

	%* Write out scoring code;
	data _Null_;
		file "&path_output.CE5_Scoring_equation.txt"; 
		set out.ce5_parameterestimates (keep=variable estimate) end=eof;
		if variable = 'Intercept' then do;
			put " " '0d'x;
			put "*****  Model  *****;" '0d'x;
			put " " '0d'x;
			%if &binary_dv = Y %then %do;
				put "Logit = 1 * ("  estimate +(-1) ") +" '0d'x;
			%end;
			%else %do;
				put "Score = 1 * ("  estimate +(-1) ") +" '0d'x;
			%end;
		end;
		else do;
			if eof then put variable "* ("  estimate +(-1) ") ;" '0d'x;
				else put variable "* ("  estimate +(-1) ") +" '0d'x;
		end;
		if eof then do;
			put " " '0d'x;
			%if &binary_dv = Y %then %do;
				put "Score=1/(1+exp(-Logit));" '0d'x;
			%end;
		end;
	run;

	proc datasets library=work nolist;
		delete mod val;
	run;
	quit;

%mend CE_Model_Lift;

**********************************************************************************************************************;