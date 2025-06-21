**********************************************************************************************************************;
*** 2.EDA,profiling and Recode ***;
**********************************************************************************************************************;

** Common ordinal and continuous processing macro **;
%macro pnum(insdn,var,nmiss,typ);
	%* Get Non-missing key statistics;
	proc means data=&insdn mean median min p1 p25 p75 p99 max NOPRINT;
		var &var;
		output out=EDA mean=var_mean median=var_median min=var_min p1=var_p1
			p25=var_p25 p75=var_p75 p99=var_p99 max=var_max / noinherit;
	run;
	%let skip=1;

	data EDA;
		set EDA;
		%* calculate upper and lower bounds;
		iqr=Max((var_p75-var_p25),(var_p99-var_p75),(var_p25-var_p1));
		var_lb=Min(Max((var_p25- 1.5*iqr),var_min),var_p1);
		var_ub=Max(Min((var_p75+ 1.5*iqr),var_max),var_p99);
		if var_lb=var_ub then do;
			var_lb=var_min;
			var_ub=var_max;
		end;
		var_mid=(var_max - var_min)/2;

		CALL SYMPUTX('var_LB' ,var_lb);
		CALL SYMPUTX('var_UB' ,var_ub);
		CALL SYMPUTX('var_median' ,var_median);

		if upcase("&&impmethod&typ") in ('MEAN','STD') then CALL
			SYMPUTX('var_miss' ,var_mean);
		else if upcase("&&impmethod&typ") in ('MEDIAN','IQR','MAD') then CALL
			SYMPUTX('var_miss' ,var_median);
		else if upcase("&&impmethod&typ") in ('RANGE') then CALL
			SYMPUTX('var_miss' ,var_min);
		else if upcase("&&impmethod&typ") in ('MIDRANGE') then CALL
			SYMPUTX('var_miss' ,var_mid);
		else if upcase("&&impmethod&typ") in ('SUM','EUCLEN','USTD','MAXABS')
			then CALL SYMPUTX('var_miss' ,0);
		else CALL SYMPUTX('skip',0);

	run;

	%* Cap_FLR before transformation;
	data NONMISSING;
		set &insdn (where=(missing(&var)=0));
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

	%let ck=%sysevalf(&nmiss - &min_size < 0);
	%if &skip=0 and %upcase(&&impmethod&typ) ^= ER %then %do;
		proc stdize data=&insdn out=dummy outstat=tmp method=&&impmethod&typ;
			var &var;
		run;

		data _NULL_;
			set tmp (where=(_type_="LOCATION"));
			call symputx("var_miss",&var);
		run;

		proc datasets nolist;
			delete tmp dummy;
		quit;
		run;
	%end;
	%else %if %upcase(&&impmethod&typ)=ER and &ck=0 %then %do;
		%* Get missing data average target;
		proc means data=&insdn mean noprint;
			var &dep_var;
			where missing(&var);
			output out=M_RR mean=Missing_RR;
		run;

		data _NULL_;
			set M_RR;
			%if %upcase(&Binary_dv)=Y %then %do;
				if Missing_RR=0 then Missing_RR=0.0001;
				else if Missing_RR=1 then Missing_RR=0.9999;
			%end;
			call symputx('Missing_RR' ,Missing_RR);
		run;

		proc datasets nolist;
			delete M_RR;
		quit;
		run;
	%end;

	%if %upcase(&Binary_dv)=Y %then %do;
		%* Run univariate logistic regression;
		ods listing close;
		ods output parameterestimates=parm association=assoc;

		proc logistic data=NONMISSING desc namelen=32;
			model &dep_var=&mod_list /PARMLABEL selection=forward STOP=1
				slentry=1;
		run;
		ods listing;
		%if %sysfunc(exist(assoc)) %then %do;
			data _NULL_;
				set assoc(keep=cvalue1 cvalue2);
				if _n_=1 then call symputx("Concordant",cvalue1);
				else if _N_=4 then call symputx("CValue",cvalue2);
			run;

			proc datasets nolist;
				delete assoc;
			quit;
			run;
		%end;
		%else %do;
			%let Concordant=. ;
			%let CValue=. ;
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

			data _NULL_;
				set SelectionSummary(keep=ModelRsquare obs=1);
				call symputx("RSquare",ModelRsquare);
			run;

		proc datasets nolist;
			delete SelectionSummary;
		quit;
		run;
	%end;
	%else %do;
		%let RSquare=. ;
	%end;
	%end;

	data _NULL_;
		set parm(obs=1);
		call symputx('Intercept',Estimate);
	run;

	data _NULL_;
		set parm (firstobs=2 obs=2);
		%if %upcase(&Binary_dv)=Y %then %do;
			call symputx('prob',ProbChiSq);
		%end;
		%else %do;
			call symputx('prob',ProbF);
		%end;
	run;
	%let ckP=%sysevalf(&prob - .05 > 0);

	data parm;
		set parm (firstobs=2 obs=2);
		length relationship $ 6 ck 3.;

		_Trans=substr(variable,1,3);
		if _Trans not in ('SQ_','SR_','LN_','IV_','EP_') or (variable="&var")
			then _Trans='';

		%if %upcase(&&impmethod&typ)=ER %then %do;
			%if %upcase(&Binary_dv)=Y %then %do;
				%*missing imputation based on missing target rate;
				%if &ck=1 or &ckP=1 %then %do;
					miss_impute=&var_median;
				%end;
				%else %do;
					if _Trans='SQ_' then
						miss_impute=SQRT(Max((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate,0));
					else if _Trans='SR_' then
						miss_impute=((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate)**2;
					else if _Trans='LN_' then
						miss_impute=exp((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate);
					else if _Trans='IV_' then
						miss_impute=-1/((log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate);
					else if _Trans='EP_' then
						miss_impute=-log(Max(-(log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate,0.00001));
					else
						miss_impute=(log(&Missing_RR/(1-&Missing_RR))-&Intercept)/Estimate;
				%end;
			%end;
			%else %do;
				%*missing imputation based on missing target rate;
				%if &ck=1 or &ckP=1 %then %do;
					miss_impute=&var_median;
				%end;
				%else %do;
					if _Trans='SQ_' then miss_impute=
						SQRT(Max((&Missing_RR-&Intercept)/Estimate,0));
					else if _Trans='SR_' then miss_impute=
						((&Missing_RR-&Intercept)/Estimate)**2;
					else if _Trans='LN_' then miss_impute=
						exp((&Missing_RR-&Intercept)/Estimate);
					else if _Trans='IV_' then miss_impute=
						-1/((&Missing_RR-&Intercept)/Estimate);
					else if _Trans='EP_' then miss_impute=
						-log(Max(-(&Missing_RR-&Intercept)/Estimate,0.00001));
					else miss_impute=(&Missing_RR-&Intercept)/Estimate;
				%end;
			%end;
			if miss_impute<&var_LB then miss_impute=&var_LB;
			else if miss_impute>&var_UB then miss_impute=&var_UB;
			call symputx("var_miss",miss_impute);
		%end;

		if Estimate>0 then sign='(+)';
		else sign='(-)';
		relationship=compress(_Trans||sign,'');

		** Output fields to append to the vars table;
		call symputx("New_var",variable);
		call symputx("Sign",sign);
		call symputx("Relation",relationship);
		%if %upcase(&Binary_dv)=Y %then %do;
			call symputx("Prob",ProbChiSq);
			call symputx("Chisq",WaldChiSq);
		%end;
		%else %do;
			call symputx("Prob",ProbF);
			call symputx("FValue",FValue);
		%end;
	run;

	%** Update vars table **;
	data vars2;
		set vars2;
		if name="&var" then do;
			var_lb=&var_lb;
			var_ub=&var_ub;
			var_median=&var_median;
			miss_impute=&var_miss;
			%if %upcase(&Binary_dv)=Y %then %do;
				Concordant=&Concordant;
				CValue=&CValue;
				Chisq=&Chisq;
				PValue=&Prob;
			%end;
			%else %do;
				RSquare=&RSquare;
				FValue=&FValue;
				PValue=&Prob;
			%end;
			if PValue <= &PValue then do;
				new_var="&Prefix"||"&New_var";
				Sign="&Sign";
				Relationship="&Relation";
			end;
		end;
	run;

	%if %sysevalf(&Pvalue-&prob>=0) %then %do;
		%if %upcase("&&stdmethod&typ") ^= "NO" %then %do;
			%* Standardization;
			data temp;
				set &insdn (keep=&var);
				_Trans=substr("&Relation",1,3);
				&New_var=&var;
				if missing(&New_var) then &New_var=&var_miss;
				%if %upcase(&&cap_flr&typ)=Y %then %do;
					&New_var=MIN(MAX(&New_var, &Var_LB), &var_UB);
				%end;
				*standard transformation;
				if _Trans='SQ_' then &New_var=&New_var **2;
				else if _Trans='SR_' then &New_var=SQRT(MAX(&New_var,0));
				else if _Trans='LN_' then &New_var=LOG(MAX(&New_var,0.00001));
				else if _Trans='IV_' then &New_var=-1/(MAX(&New_var,0.00001));
				else if _Trans='EP_' then &New_var=-EXP(MIN(-&New_var,0));
			run;

			proc stdize data=temp out=dummy outstat=std_tmp method=
				&&stdmethod&typ;
				var &New_var;
			run;

			data _NULL_;
				set std_tmp;
				if _type_='LOCATION' then call symputx("loc",&new_var);
				if _type_='SCALE' then call symputx("scale",&new_var);
			run;

			%** Update vars table **;
			data vars2;
				set vars2;
				if name="&var" then do;
					Location=&loc;
					Scale=&scale;
				end;
			run;

			proc datasets nolist;
				delete temp dummy std_tmp;
			quit;
			run;
		%end;

		data _NULL_;
			set vars2 (where=(name="&var"));
			_Trans=substr(relationship,1,3);
			if missing(label) then lab2=name;
			else lab2=label;

			%if &typ=O %then %do;
				FILE "&Path_output.CE2_Ordinal_Var_Recode.txt" mod;
			%end;
			%else %do;
				FILE "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
			%end;
			PUT ' ' '0d'x;
			PUT "*** RECODE " name ": " label " ***;" '0d'x;
			* Missing imputation;
			Put "IF missing(" name +(-1) ") THEN  &PREFIX" name " = "
				Miss_Impute " ; " '0d'x;
			* Valid values;
			PUT @2 "ELSE &PREFIX" name " = " name ";" '0d'x;
			* Capping and flooring;
			%if %upcase(&&cap_flr&typ)=Y %then %do;
				Put "&PREFIX" name " = MIN(MAX(&PREFIX" name ", " Var_LB +(-1)
					"), " var_UB +(-1) ");" '0d'x;
			%end;
			* Untransformed variable label;
			%if &typ=O %then %do;
				PUT @2 "Label &PREFIX" name ' = "' lab2 ': Ordinal Recode ' Sign
					'";' '0d'x;
			%end;
			%else %do;
				PUT @2 "Label &PREFIX" name ' = "' lab2 ': Continuous Recode '
					Sign '";' '0d'x;
			%end;

			*standard transformation;
			if _Trans='SQ_' then do;
				Put new_var "= &PREFIX" name "**2;" '0d'x;
				Put @2 " Label " new_var '= "' lab2 'SQUARE ' Sign '";' '0d'x;
			end;
			else if _Trans='SR_' then do;
				Put new_var " = SQRT(MAX(&PREFIX" name ",0));" '0d'x;
				Put @2 " Label " new_var '= "' lab2 'SQRT ' Sign '";' '0d'x;
			end;
			else if _Trans='LN_' then do;
				Put new_var "= LOG(MAX(&PREFIX" name ",0.00001));" '0d'x;
				Put @2 " Label " new_var '= "' lab2 'LOG ' Sign '";' '0d'x;
			end;
			else if _Trans='IV_' then do;
				Put new_var "= -1/(MAX(&PREFIX" name ",0.00001));" '0d'x;
				Put @2 " Label " new_var '= "' lab2 'NEGATIVE INVERSE ' Sign
					'";' '0d'x;
			end;
			else if _Trans='EP_' then do;
				Put new_var "= -EXP(MIN(-&PREFIX" name ",0));" '0d'x;
				Put @2 " Label " new_var '= "' lab2 'EXPONENTIAL ' Sign '";'
					'0d'x;
			end;
			* Standardization;
			%if %upcase("&&stdmethod&typ") ^= "NO" %then %do;
				if (location ^= . and scale ^= .) then put new_var "= (" new_var
					"- (" location +(-1) ") ) / " scale ";" '0d'x;
			%end;

		run;
	%end;

	proc datasets nolist;
		delete parm nonmissing;
	quit;
	run;
%mend pnum;

** Ordinal and continuous profiling macro 1 **;
%macro prof1(insdn,var);
	proc summary data=&insdn nway missing;
		var &dep_var;
		class &var;
		output out=prof (drop=_type_ rename=_freq_=xcount) mean=xmean;
	run;

	data prof;
		set prof;
		length xcategory $256.;
		if missing(&var) then xcategory='Missing';
		else xcategory=trim(left(put(&var,best8.)));
	run;
%mend;

** Ordinal and continuous profiling macro 2 **;
%macro prof2(insdn,var);
	%if %upcase(&equal_dist)=Y %then %do;
		data _NULL_;
			set EDA;
			range=(var_p99 - var_p1)/&num_category;
			call symputx('cut_lo',var_p1);
			call symputx('cut_hi',var_p99);
			call symputx('range',range);
		run;
		%put &range;

		data tmp;
			set &insdn (keep=&dep_var &var);
			if missing(&var) then bin=.;
			else if &var < &cut_lo then bin=1;
			else if &var >= &cut_hi then bin=input("&num_category",best10.);
			else do;
				do k=1 to &num_category;
					if &var>= &cut_lo+(k-1)*&range and &var< &cut_lo+k*&range
						then bin=k;
				end;
			end;
		run;
	%end;
	%else %do;
		proc rank data=&insdn (keep=&dep_var &var) out=tmp ties=High
			group=&num_category;
			var &var;
			ranks bin;
		run;
	%end;

	proc summary data=tmp nway missing;
		var &var &dep_var;
		class bin;
		output out=prof (drop=_type_ rename=_freq_=xcount) min(&var)=lo
			max(&var)=hi mean(&dep_var)=xmean;
	run;

	data prof (drop=tag);
		set prof end=eof;
		length xcategory $256.;
		retain tag 0;
		if missing(bin) then do;
			xcategory='Missing';
			tag=1;
		end;
		else if _N_=1 or (_N_=2 and tag=1) then xcategory="Low to " ||
			trim(left(put(hi,best8.)));
		else if eof then xcategory=trim(left(put(lo,best8.))) || " to High";
		else xcategory=trim(left(put(lo,best8.))) || " to " ||
			trim(left(put(hi,best8.)));
	run;

	proc datasets nolist;
		delete tmp eda;
	run;
%mend;

** Ordinal and continuous profiling macro 3 **;
%macro prof3(typ);
	proc sql;
		create table prof2 as select a.name as variable, a.label, b.xcategory,
			b.xcount, b.xmean from vars2 a, prof b where a.name="&var";
	quit;

	data prof3 (drop=xcount xmean xcategory);
		set prof2 end=eof;
		length category $256.;
		length type $10.;
		length star $8.;
		if "&typ"="C" then type='Continuous';
		else type='Ordinal';
		if _N_=1 then do;
			category="Overall";
			Average_DV=&overall_avg;
			Count=&nobs;
			Percent=1;
			index=100;
			output;
		end;
		category=xcategory;
		count=xcount;
		percent=xcount / &nobs;
		Average_DV=xmean;
		index=(Average_DV / &overall_avg)*100;
		if index >= 110 then star='* (+)';
		else if index > 100 then star='  (+)';
		else if index <= 90 then star='* (-)';
		else if index <= 100 then star='  (-)';
		else star='  (0)';
		output;
	run;

	data out.CE2_profile;
		set out.CE2_profile prof3;
	run;

	proc datasets nolist;
		delete prof prof2 prof3;
	run;
%mend;

** Binary variable processing macro **;
%macro pbin(insdn,var,typ);
	%** Get missing count **;
	proc sql noprint;
		select max(sum(missing(&var)),0) into : miss_ck from &insdn;
	quit;
	%put processing binary variable &var with type &typ;
	%if %sysevalf((&nobs*&missrate)>=&miss_ck) %then %put missing count is
		&miss_ck;
	%else %put missing count is &miss_ck which is too high;
	%** If missing count is acceptable **;
	%if %sysevalf((&nobs*&missrate)>=&miss_ck) %then %do;
		%** Get counts of values that should be coded as 1 **;
		%if &typ=1 %then %do; ** Numeric variables **;
			proc sql noprint;
				select sum(case when &var=1 then 1 else 0 end) into : cnt_ck
					from &insdn;
			quit;
		%end;
		%else %do; %** Character variables **;
			proc sql noprint;
				select sum(case when strip(&var) in ('1','Y','y') then 1 else 0
					end) into : cnt_ck from &insdn;
			quit;
		%end;
		%if %sysevalf((&nobs*(1-&concrate)) <= &cnt_ck and &cnt_ck <=
			(&nobs*&concrate)) %then %do;
			%** If count is okay, write out code to file **;
			data _null_;
				set vars (where=(name="&var"));
				FILE "&Path_output.CE2_Binary_Var_Recode.txt" mod;
				PUT ' ' '0d'x;
				PUT "*** RECODE " name ": " label " ***;" '0d'x;
				%if &typ=1 %then %do;
					PUT "if " name " = 1 then &PREFIX" name " = 1; else &PREFIX"
						name " = 0;" '0d'x;
				%end;
				%else %do;
					PUT "if " name " in ('1','Y','y') then &PREFIX" name
						" = 1; else &PREFIX" name " = 0;" '0d'x;
				%end;
				if missing(label) then PUT @2 "Label &PREFIX" name ' = "' name
					': Binary Recode";' '0d'x;
				else PUT @2 "Label &PREFIX" name ' = "' label
					': Binary Recode";' '0d'x;
			run;
			%if %upcase(&profiling)=Y %then %do;
				%if &typ=1 %then %do;
					data tmp1;
						set &insdn (keep=&dep_var &var);
						if &var=1 then newvar=1;
						else newvar=0;
					run;
				%end;
				%else %do;
					data tmp1;
						set &insdn;
						if &var in ('1','Y','y') then newvar=1;
						else newvar=0;
					run;
				%end;

				proc sql;
					create table tmp2 as select a.name as variable, a.label,
						b.newvar, count(*) as xcount, mean(b.&dep_var) as xmean
						from vars a, tmp1 b where a.name="&var" group by a.name,
						a.label, b.newvar order by a.name, a.label, b.newvar;
				quit;

				data tmp3 (drop=xcount xmean newvar);
					set tmp2 end=eof;
					by newvar;
					length type $10.;
					length category $256.;
					length star $8.;
					type='Binary';
					if _N_=1 then do;
						category="Overall";
						Average_DV=&overall_avg;
						Count=&nobs;
						Percent=1;
						index=100;
						output;
					end;
					count=xcount;
					percent=xcount / &nobs;
					Average_DV=xmean;
					index=(Average_DV / &overall_avg)*100;
					%if &typ=1 %then %do;
						if newvar=0 then category="Missing,0";
						else category="1";
					%end;
					%else %do;
						if newvar=0 then category="Missing,0,N";
						else category="1,Y";
					%end;
					if index >= 110 then star='* (+)';
					else if index > 100 then star='  (+)';
					else if index <= 90 then star='* (-)';
					else if index <= 100 then star='  (-)';
					else star='  (0)';
					output;
				run;

				data out.CE2_profile;
					set out.CE2_profile tmp3;
				run;

				proc datasets nolist;
					delete tmp: ;
					run;
					%end;
					%end;
					%else %do;
						%** If variable is too concentrated, write note to log **;
						%put variable &var is too concentrated to use;
					%end;

				%** Add counts to summary file **;
				data vars;
					set vars;
					if name="&var" then do;
						miss_cnt=&miss_ck;
						ones_cnt=&cnt_ck;
						%if %sysevalf((&nobs*(1-&concrate)) <= &cnt_ck and
							&cnt_ck <= (&nobs*&concrate)) %then %do;
							good=1;
						%end;
					end;
				run;
			%end;
			%** If missing count is too high **;
			%else %do;
				data vars;
					set vars;
					if name="&var" then do;
						miss_cnt=&miss_ck;
					end;
				run;
			%end;
			%mend pbin;

		** Binary processing control macro **;
		%macro bin_cntl(insdn);
			%** Initialize code file **;
			data _null_;
				FILE "&Path_output.CE2_Binary_Var_Recode.txt" LRECL=256;
				PUT
					"*********************************************************************************;"
					'0d'x/
					"****                        CE2_BINARY_VAR_RECODE.TXT                        ****;"
					'0d'x/
					"*********************************************************************************;"
					'0d'x;
			run;

			%** Get variables to process **;
			proc contents data=&insdn (keep=&binvar) out=vars (keep=varnum name
				type length label) noprint;
			run;

			proc sql noprint;
				select count(*) into : ck from vars;
			quit;
			%put number of binary variables=&ck;

			%** If there are binary variables to process **;
			%if &ck > 0 %then %do;
				%** Add new variables to summary file **;
				data vars;
					set vars;
					length miss_cnt ones_cnt good 8.;
				run;

				%** Create macro variables to control processing **;
				data _null_;
					set vars (where=(length(strip(name)) +
						length(strip("&prefix")) <= 32)) end=eof;
					if eof then call symputx("VARCNT",_N_);
					call symputx("IV"|| trim(left(put(_N_,4.))) ,name);
					call symputx("typ"|| trim(left(put(_N_,4.))) ,type);
				run;
				%** Loop through for each variable **;
				%do _I_=1 %to &varcnt;
					%pbin(&insdn,&&iv&_I_,&&typ&_I_);
				%end;

				%** Create new variable name for surviving variables **;
				data out.CE2_binary_vars (drop=good);
					set vars;
					length new_var $32.;
					if good=1 then new_var="&prefix"||name;
				run;
				%** Write out list of variables to keep to code file **;
				%let rck=0;

				data _NULL_;
					set out.CE2_binary_vars (where=(missing(new_var)=0))
						end=eof;
					m=mod(_N_,7);
					length v $256.;
					retain v ;
					if m=1 then v=new_var;
					else v=strip(v) || " " || strip(new_var);
					FILE "&Path_output.CE2_Binary_Var_Recode.txt" mod;
					if _N_=1 then do;
						PUT ' ' '0d'x;
						PUT ' ' '0d'x;
						PUT '%LET KEEPLIST_B =';
					end;
					if m=0 or eof then PUT v '0d'x;
					if eof then PUT ';';
					if eof then call symputx("rck" ,_N_);
				run;
				%put number of recoded binary variables=&rck;

				%* Capture old and new variable name;
				data tmp (keep=orig_var orig_label variable);
					length variable orig_var $32.;
					length orig_label $256.;
					set out.CE2_binary_vars (where=(missing(variable)=0)
						rename=(name=orig_var label=orig_label
						new_var=variable));
				run;

				data keep_vars;
					set keep_vars tmp;
				run;

				proc datasets nolist;
					delete vars tmp;
				quit;
				run;
			%end;
			%** If there are no binary variables to process **;
			%else %do;
				data _null_;
					FILE "&Path_output.CE2_Binary_Var_Recode.txt" mod;
					PUT ' ' '0d'x;
					PUT '%LET KEEPLIST_B = ;' '0d'x;
				run;
			%end;
		%mend bin_cntl;

	** Nominal variable processing macro **;
	%macro pnom(insdn,var,typ);
		%** Summarize **;
		proc sql;
			create table tmp as select &var, count(*) as dcount, mean(&dep_var)
				as dmean, var(&dep_var) as dvar from &insdn group by &var order
				by dmean;
		quit;

		proc sql noprint;
			select sum(case when missing(&var) then dcount else 0 end) into :
				misscnt from tmp;
			select max(dcount) into : maxcnt from tmp;
			select count(*) into : unqcnt from tmp;
		quit;
		%** Test **;
		%put processing nominal variable &var with type &typ;
		%if %sysevalf((&nobs*&missrate)<&misscnt) %then %put missing count is
			&misscnt which is too high;
		%else %do;
			%put missing count is &misscnt;
			%if %sysevalf((&nobs*&concrate)<&maxcnt) %then %put variable &var is
				too concentrated to use;
			%else %do;
				%if %sysevalf(&valcnt^=0 and &valcnt<=&unqcnt) %then %put
					variable &var has too many values to use;
				%else %do; %** Variable is usable;
					%** Collapse values based on counts;
					DATA tmp1;
						SET tmp END=eof;
						LENGTH tcount tmean tvar tgroup cumcnt tgrpcnt 8.;
						RETAIN tcount tmean tvar tgroup cumcnt tgrpcnt;
						cumcnt + dcount;
						IF _n_=1 THEN DO;
							tcount=dcount;
							tmean=dmean;
							tvar=dvar;
							tgroup=1;
							tgrpcnt=1;
						END;
						ELSE DO;
							IF tcount <= &minbinn or (cumcnt - dcount) >= (&nobs
								- &minbinn) THEN DO;
								tmean=((tcount * tmean) + (dcount *
									dmean))/(tcount + dcount);
								tvar=(((tcount - 1) * tvar) + ((dcount - 1) *
									dvar))/(tcount + dcount - 2);
								tcount=tcount + dcount;
								tgroup=tgroup;
								tgrpcnt=tgrpcnt + 1;
							END;
							ELSE DO;
								tcount=dcount;
								tmean=dmean;
								tvar=dvar;
								tgroup=tgroup+1;
								tgrpcnt=1;
							END;
						END;
						KEEP &var dcount dmean dvar tgroup tcount tmean tvar
							tgrpcnt;
					RUN;

					%** Create dataset with final record for group only;
					data tmp2;
						set tmp1 (drop=&var dcount dmean dvar);
						by tgroup;
						if last.tgroup;
					run;

					%** Calculate bonferroni adjustment on alpha value;
					%IF %UPCASE(&bonfer)=Y %THEN %DO;
						PROC SQL NOPRINT;
							SELECT COUNT(DISTINCT group) INTO :ncomps FROM tmp2;
						QUIT;
						%LET ncomps=%SYSEVALF(&ncomps - 1);
						%IF &ncomps > 1 %THEN %LET ftalpha=
							%SYSEVALF(&talpha./&ncomps);
						%ELSE %LET ftalpha=&talpha;
					%END;
					%ELSE %LET ftalpha=&talpha;

					%** Collapse values based on variance;
					DATA tmp3;
						SET tmp2 END=eof;
						LENGTH row fcount fmean fvar fgroup grpcnt 8.;
						RETAIN fcount fmean fvar fgroup grpcnt;
						row=_N_;
						IF _n_=1 THEN DO;
							fcount=tcount;
							fmean=tmean;
							fvar=tvar;
							fgroup=1;
							grpcnt=tgrpcnt;
						END;
						ELSE DO;
							pvar=(((fcount - 1) * fvar) + ((tcount - 1) *
								tvar))/(fcount + tcount - 2);
							t=(fmean - tmean)/SQRT(pvar * ((1/fcount) +
								(1/tcount)));
							df=(fcount + tcount - 2);
							prob=(1-PROBT(abs(t),df));
							IF prob <= &ftalpha THEN DO;
								fcount=tcount;
								fmean=tmean;
								fvar=tvar;
								fgroup=fgroup+1;
								grpcnt=tgrpcnt;
							END;
							ELSE DO;
								fmean=((fcount * fmean) + (tcount *
									tmean))/(fcount + tcount);
								fvar=(((fcount - 1) * fvar) + ((tcount - 1) *
									tvar))/(fcount + tcount - 2);
								fcount=fcount + tcount;
								fgroup=fgroup;
								grpcnt=grpcnt + tgrpcnt;
							END;
						END;
						KEEP row tgroup fgroup fcount fmean grpcnt;
					RUN;

					PROC SORT DATA=tmp3 out=tmp4;
						BY fgroup DESCENDING row;
					RUN;

					DATA tmp4 (drop=fcount fmean row grpcnt);
						SET tmp4;
						LENGTH xmean xcount group_size 8. name $32.;
						RETAIN xmean xcount group_size;
						BY fgroup DESCENDING row;
						IF FIRST.fgroup THEN DO;
							xmean=fmean;
							xcount=fcount;
							group_size=grpcnt;
						END;
						Index=(xmean/&overall_Avg)*100 ;
						diff=abs(xmean-&overall_Avg);
						name="&var";
					RUN;

					%* Sort to keep smaller groups first;
					PROC SORT DATA=tmp4;
						BY group_size DESCENDING diff xcount ;
					RUN;

					DATA _NULL_;
						SET tmp4 END=eof;
						IF eof THEN call symputx('last_group',fgroup);
					RUN;

					%* Combine data;
					proc sql;
						create table tmp5 as select b.name, c.label, a.&var,
							b.fgroup, a.dmean, a.dcount, b.xmean, b.xcount,
							b.diff, b.Index, b.group_size from tmp1 a, tmp4 b,
							vars c where a.tgroup=b.tgroup and b.name=c.name
							order by b.group_size, b.diff desc, b.xcount;
					quit;

					%* Write out recode code;
					%if %upcase(&nom_method)=BINARY %then %do;
						DATA _null_;
							SET tmp5 (where=(fgroup ^= &last_group)) END=eof;
							BY group_size DESCENDING diff ;
							FILE "&Path_output.CE2_Nominal_Var_Recode.txt" mod;

							IF _n_=1 THEN DO;
								PUT ' ' '0d'x;
								PUT "*** RECODE " name ": " label "***;" '0d'x;
							END;
							%if &typ=1 %then %do;
								IF FIRST.diff THEN PUT "&PREFIX" name +(-1) "_X"
									fgroup "= " name "in (" &var +(-1)@;
								IF FIRST.diff NE 1 THEN PUT "," &var +(-1)@;
								IF LAST.diff THEN PUT ");" '0d'x;
							%end;
							%else %do;
								IF FIRST.diff THEN PUT "&PREFIX" name +(-1) "_X"
									fgroup "= " name "in ('" &var +(-1)@;
								IF FIRST.diff NE 1 THEN PUT "','" &var +(-1)@;
								IF LAST.diff THEN PUT "');" '0d'x;
							%end;
						RUN;

						DATA _null_;
							SET tmp5 (where=(fgroup ^= &last_group)) END=eof;
							BY group_size DESCENDING diff ;
							FILE "&Path_output.CE2_Nominal_Var_Recode.txt" mod;

							IF FIRST.diff THEN DO;
								IF missing(label) then PUT @2 "Label &PREFIX"
									name +(-1) "_X" fgroup '= "' name
									': Values ' &var +(-1)@;
								ELSE PUT @2 "Label &PREFIX" name +(-1) "_X"
									fgroup '= "' label ': Values ' &var +(-1)@;
							END;
							IF FIRST.diff NE 1 THEN PUT "," &var +(-1)@;
							IF LAST.diff THEN PUT '";' '0d'x;
						RUN;

						data nv (keep=name new_var);
							set tmp5 (where=(fgroup ^= &last_group));
							by group_size DESCENDING diff ;
							if first.diff;
							length new_var $32.;
							new_var="&prefix" || strip(name) || "_X" ||
								strip(left(put(fgroup,best2.)));
						run;

						data new_vars;
							set new_vars nv;
						run;
					%end;
					%else %do;
						%if %upcase(&nom_method)=INDEX %then %let val=index;
						%else %let val=xmean;

						proc sort data=tmp5;
							by fgroup group_size &var;
						run;

						proc sql noprint;
							select max(fgroup) into : last_group from tmp5;
						quit;

						DATA _null_;
							SET tmp5 END=eof;
							BY fgroup group_size;
							FILE "&Path_output.CE2_Nominal_Var_Recode.txt" mod;

							IF _n_=1 THEN DO;
								IF fgroup=&last_group THEN STOP;
								ELSE DO;
									PUT ' ' '0d'x;
									PUT "*** RECODE " name ": " label "***;"
										'0d'x;
								END;
							END;

							IF fgroup=&last_group THEN DO;
								IF first.group_size THEN PUT "else &PREFIX" name
									"= " &val ";" '0d'x;
							END;
							ELSE DO;
								%if &typ=1 %then %do;
									IF _n_=1 AND FIRST.group_size THEN PUT "if "
										name "in (" &var +(-1)@;
									ELSE IF FIRST.group_size THEN PUT "else if "
										name "in (" &var +(-1)@;
									IF FIRST.group_size NE 1 THEN PUT "," &var
										+(-1)@;
									IF LAST.group_size THEN PUT ") then &PREFIX"
										name "= " &val ";" '0d'x;
								%end;
								%else %do;
									IF _n_=1 AND FIRST.group_size THEN PUT "if "
										name "in ('" &var +(-1)@;
									ELSE IF FIRST.group_size THEN PUT "else if "
										name "in ('" &var +(-1)@;
									IF FIRST.group_size NE 1 THEN PUT "','" &var
										+(-1)@;
									IF LAST.group_size THEN PUT
										"') then &PREFIX" name "= " &val ";"
										'0d'x;
								%end;
							END;
						RUN;

						DATA _null_;
							SET tmp5 (where=(fgroup ^= &last_group) obs=1)
								END=eof;
							FILE "&Path_output.CE2_Nominal_Var_Recode.txt" mod;
							IF missing(label) then PUT @2 "Label &PREFIX" name
								' = "' name ': Nominal Recode";' '0d'x;
							ELSE PUT @2 "Label &PREFIX" name ' = "' label
								': Nominal Recode";' '0d'x;
						RUN;

						data nv (keep=name new_var);
							set tmp5 (where=(fgroup ^= &last_group) obs=1);
							length new_var $32.;
							new_var="&prefix" || strip(name);
						run;

						data new_vars;
							set new_vars nv;
						run;
					%end;
					%if %upcase(&Profiling)=Y %then %do;
						proc sort data=tmp5;
							by fgroup group_size &var;
						run;

						data tmp6 (keep=type Variable label category count
							percent Average_DV index star);
							set tmp5 (rename=(name=Variable)) end=eof;
							by fgroup;
							length type $10.;
							length category $256.;
							length star $8.;
							retain category;
							type='Nominal';
							if _N_=1 then do;
								category="Overall";
								Average_DV=&overall_avg;
								Count=&nobs;
								Percent=1;
								index=100;
								output;
							end;
							count=xcount;
							percent=xcount / &nobs;
							Average_DV=xmean;
							Index=(Average_DV / &overall_avg)*100;
							if first.fgroup and missing(&var) then category=
								"Missing";
							else if first.fgroup then category=strip(&var);
							else category=strip(category) || "," || strip(&var);
							if index >= 110 then star='* (+)';
							else if index > 100 then star='  (+)';
							else if index <= 90 then star='* (-)';
							else if index <= 100 then star='  (-)';
							else star='  (0)';
							if last.fgroup then output;
						run;

						data out.CE2_profile;
							set out.CE2_profile tmp6;
						run;
					%end;
				%end;
			%end;
		%end;

		%** Add counts to summary file **;
		data vars;
			set vars;
			if name="&var" then do;
				unique_cnt=&unqcnt;
				miss_cnt=&misscnt;
				max_cnt=&maxcnt;
			end;
		run;

		proc datasets nolist;
			delete tmp: nv;
		quit;
		run;

	%mend pnom;

** Nominal processing control macro **;
%macro nom_cntl(insdn);
	%** Initialize code file **;
	data _null_;
		FILE "&Path_output.CE2_Nominal_Var_Recode.txt" LRECL=256;
		PUT
			"*********************************************************************************;"
			'0d'x/
			"****                        CE2_NOMINAL_VAR_RECODE.TXT                       ****;"
			'0d'x/
			"*********************************************************************************;"
			'0d'x;
	run;

	%** Get variables to process **;
	proc contents data=&insdn (keep=&nomvar) out=vars (keep=varnum name type
		length label) noprint;
	run;

	proc sql noprint;
		select count(*) into : ck from vars;
	quit;
	%put number of nominal variables=&ck;

	%** If there are nominal variables to process **;
	%if &ck > 0 %then %do;
		%** Add new variables to summary file **;
		data vars;
			set vars;
			length unique_cnt miss_cnt max_cnt 8.;
		run;

		%if %upcase(&nom_method)=BINARY %then %let maxlen=29;
		%else %let maxlen=32;

		%** Create macro variables to control processing **;
		data _null_;
			set vars (where=(length(strip(name)) + length(strip("&prefix")) <=
				&maxlen)) end=eof;
			if eof then call symputx("VARCNT",_N_);
			call symputx("IV"|| trim(left(put(_N_,4.))) ,name);
			call symputx("typ"|| trim(left(put(_N_,4.))) ,type);
		run;

		data new_vars;
			set _NULL_;
			length name new_var $32.;
		run;
		%** Loop through for each variable **;
		%do _I_=1 %to &varcnt;
			%pnom(&insdn,&&iv&_I_,&&typ&_I_);
		%end;

		%** Create new variable name for surviving variables **;
		proc sort data=new_vars;
			by name;
		run;

		data nv;
			set new_vars;
			by name;
			length new_vars $256.;
			retain new_vars;
			if first.name then new_vars=new_var;
			else new_vars=strip(new_vars)||","||strip(new_var);
			if last.name then output;
		run;

		proc sql;
			create table out.CE2_nominal_vars as select a.*, b.new_vars from
				vars a left join nv b on a.name=b.name;
		quit;
		%** Write out list of variables to keep to code file **;
		%let rck=0;

		data _NULL_;
			set new_vars (where=(missing(new_var)=0)) end=eof;
			m=mod(_N_,7);
			length v $256.;
			retain v ;
			if m=1 then v=new_var;
			else v=strip(v) || " " || strip(new_var);
			FILE "&Path_output.CE2_Nominal_Var_Recode.txt" mod;
			if _N_=1 then do;
				PUT ' ' '0d'x;
				PUT ' ' '0d'x;
				PUT '%LET KEEPLIST_N =';
			end;
			if m=0 or eof then PUT v '0d'x;
			if eof then PUT ';';
			if eof then call symputx("rck" ,_N_);
		run;
		%put number of recoded nominal variables=&rck;

		%* Capture old and new variable name;
		proc sql;
			create table tmp as select a.name as orig_var, a.label as
				orig_label, b.new_var as variable from vars a, new_vars b where
				a.name=b.name;
		quit;

		data keep_vars;
			set keep_vars tmp;
		run;

		proc datasets nolist;
			delete vars new_vars nv tmp;
		quit;
		run;
	%end;
	%** If there are no nominal variables to process **;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Nominal_Var_Recode.txt" mod;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_N = ;' '0d'x;
		run;
	%end;
%mend nom_cntl;

** Ordinal variable processing macro **;
%macro pord(insdn,var,nmiss);
	%** Summarize **;
	proc sql;
		create table tmp as select &var, count(*) as count from &insdn group by
			&var;
	quit;

	proc sql noprint;
		select max(count) into : maxcnt from tmp;
		select count(*) into : unqcnt from tmp;
	quit;

	%* check if dep_var has the constant value in nonmissing part;
	proc sql noprint;
		select min(&dep_var)=max(&dep_var) into : cons_dep from &insdn where
			missing(&var)=0;
	quit;

	proc datasets nolist;
		delete tmp;
	quit;
	run;

	%* add to vars file;
	data vars2;
		set vars2;
		if name="&var" then do;
			unique_cnt=&unqcnt;
			max_cnt=&maxcnt;
		end;
	run;

	%** Test **;
	%if %sysevalf((&nobs*&concrate)<&maxcnt) %then %put variable &var is too
		concentrated to use;
	%else %if &cons_dep=1 %then %put &var has constant dependent value in
		nonmissing part and will be excluded;
	%else %do;
		%pnum(&insdn,&var,&nmiss,O);
		%** Profiling **;
		%if %upcase(&profiling)=Y %then %do;
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
		FILE "&Path_output.CE2_Ordinal_Var_Recode.txt" LRECL=256;
		PUT
			"*********************************************************************************;"
			'0d'x/
			"****                       CE2_ORDINAL_VAR_RECODE.TXT                        ****;"
			'0d'x/
			"*********************************************************************************;"
			'0d'x;
	run;

	%** Get variables to process **;
	proc contents data=&insdn (keep=&ordvar) out=vars (keep=varnum name type
		length label) noprint;
	run;

	proc sql noprint;
		select max(type) into : ck_type from vars;
	quit;
	%if &ck_type=2 %then %do;
		proc sql noprint;
			select name into : bad separated by ',' from vars where type=2;
		quit;
		%put The following character variables were dropped from the ordinal
			variable list:;
		%put &bad;

		proc sql noprint;
			select name into : ordvar separated by ' ' from vars where type=1;
		quit;
	%end;

	data vars;
		set vars (where=(type=1));
	run;

	proc sql noprint;
		select count(*) into : ck from vars;
	quit;
	%put number of ordinal variables=&ck;

	%** If there are ordinal variables to process **;
	%if &ck > 0 %then %do;
		%** Get missing counts to check variable validity **;
		proc means data=&insdn noprint ;
			var &ordvar;
			output out=nmiss (drop=_type_ _freq_) nmiss= ;
		run;

		proc transpose data=nmiss out=nmiss name=name;
		run;

		%** Combine **;
		proc sql;
			create table vars2 as select a.*, b.col1 as miss_cnt from vars a
				left join nmiss b on a.name=b.name;
		quit;

		data vars2;
			set vars2;
			length new_var $32.;
			length sign $3.;
			length relationship $6.;
			length max_cnt unique_cnt 8.;
		run;

		%if %upcase(&transformationO)=Y %then %let maxlen=29;
		%else %let maxlen=32;

		%** Create macro variables to control processing **;
		data _null_;
			set vars2 (where=(length(strip(name)) + length(strip("&prefix")) <=
				&maxlen)) end=eof;
			if eof then call symput("VARCNT",_N_);
			call symput("IV"|| trim(left(put(_N_,4.))) ,name);
			call symput("miss"|| trim(left(put(_N_,4.))) ,miss_cnt);
		run;
		%** Loop through for each variable **;
		%do _I_=1 %to &varcnt;
			%put processing ordinal variable &&iv&_I_;
			%if %sysevalf((&nobs*&missrate)<&&miss&_I_) %then %put missing count
				is &&miss&_I_ which is too high;
			%else %do;
				%put missing count is &&miss&_I_ ;
				%pord(&insdn,&&iv&_I_,&&miss&_I_);
			%end;
		%end;

		%** Create permanent variable summary file **;
		data out.CE2_ordinal_vars;
			set vars2;
		run;
		%** Write out list of variables to keep to code file **;
		%let rck=0;

		data _NULL_;
			set out.CE2_ordinal_vars (where=(missing(new_var)=0)) end=eof;
			m=mod(_N_,7);
			length v $256.;
			retain v ;
			if m=1 then v=new_var;
			else v=strip(v) || " " || strip(new_var);
			FILE "&Path_output.CE2_Ordinal_Var_Recode.txt" mod;
			if _N_=1 then do;
				PUT ' ' '0d'x;
				PUT ' ' '0d'x;
				PUT '%LET KEEPLIST_O =';
			end;
			if m=0 or eof then PUT v '0d'x;
			if eof then PUT ';';
			if eof then call symputx("rck" ,_N_);
		run;
		%put number of recoded ordinal variables=&rck;

		%* Capture old and new variable name;
		data tmp (keep=orig_var orig_label variable);
			length variable orig_var $32.;
			length orig_label $256.;
			set out.CE2_ordinal_vars (where=(missing(variable)=0)
				rename=(name=orig_var label=orig_label new_var=variable));
		run;

		data keep_vars;
			set keep_vars tmp;
		run;

		proc datasets nolist;
			delete vars vars2 tmp;
		quit;
		run;
	%end;
	%** If there are no ordinal variables to process **;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Ordinal_Var_Recode.txt" mod;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_O = ;' '0d'x;
		run;
	%end;
%mend ord_cntl;

** Continuous variable processing macro **;
%macro pcont(insdn,var,nmiss);
	%* check if dep_var has the constant value in nonmissing part;
	proc sql noprint;
		select min(&dep_var)=max(&dep_var) into : cons_dep from &insdn where
			missing(&var)=0;
	quit;

	%if &cons_dep=1 %then %put &var has constant dependent value in nonmissing
		part and will be excluded;
	%else %do;
		%pnum(&insdn,&var,&nmiss,C);
		%** Profiling **;
		%if %upcase(&profiling)=Y %then %do;
			%prof2(&insdn,&var);
			%prof3(C);
		%end;
	%end;

%mend pcont;

** Continuous processing control macro **;
%macro cont_cntl(insdn);
	%** Initialize code file **;
	data _null_;
		FILE "&Path_output.CE2_Continuous_Var_Recode.txt" LRECL=256;
		PUT
			"*********************************************************************************;"
			'0d'x/
			"****                      CE2_CONTINUOUS_VAR_RECODE.TXT                      ****;"
			'0d'x/
			"*********************************************************************************;"
			'0d'x;
	run;

	%** Get variables to process **;
	proc contents data=&insdn (keep=&contvar) out=vars (keep=varnum name type
		length label) noprint;
	run;

	proc sql noprint;
		select max(type) into : ck_type from vars;
	quit;
	%if &ck_type=2 %then %do;
		proc sql noprint;
			select name into : bad separated by ',' from vars where type=2;
		quit;
		%put The following character variables were dropped from the continuous
			variable list:;
		%put &bad;

		proc sql noprint;
			select name into : contvar separated by ' ' from vars where type=1;
		quit;
	%end;

	data vars;
		set vars (where=(type=1));
	run;

	proc sql noprint;
		select count(*) into : ck from vars;
	quit;
	%put number of continuous variables=&ck;

	%** If there are continuous variables to process **;
	%** If there are continuous variables to process **;
	%if &ck > 0 %then %do;
		%** Get statistics to check variable validity **;
		proc means data=&insdn noprint ;
			var &contvar;
			output out=nmiss nmiss= ;
		run;

		proc transpose data=nmiss out=nmiss name=name;
		run;

		%** Combine **;
		proc sql;
			create table vars2 as select a.*, b.col1 as miss_cnt from vars a
				left join nmiss b on a.name=b.name;
		quit;

		%** Add new variables to summary file **;
		data vars2;
			set vars2;
			length P&p_lo P&p_hi 8.;
			length new_var $32.;
			length sign $3.;
			length relationship $6.;
		run;

		proc datasets nolist;
			delete nmiss;
		quit;
		run;

		%if %upcase(&transformationC)=Y %then %let maxlen=29;
		%else %let maxlen=32;

		%** Create macro variables to control processing **;
		data _null_;
			set vars2 (where=(length(strip(name)) + length(strip("&prefix")) <=
				&maxlen)) end=eof;
			if eof then call symputx("VARCNT",_N_);
			call symputx("IV"|| trim(left(put(_N_,4.))) ,name);
			call symputx("miss"|| trim(left(put(_N_,4.))) ,miss_cnt);
		run;
		%** Loop through for each variable **;
		%do _I_=1 %to &varcnt;
			%put processing continuous variable &&iv&_I_;
			%if %sysevalf((&nobs*&missrate)<&&miss&_I_) %then %put missing count
				is &&miss&_I_ which is too high;
			%else %do;
				%put missing count is &&miss&_I_ ;

				proc means data=&insdn noprint ;
					var &&iv&_I_;
					output out=st P&p_lo=P_lo P&p_hi=P_hi ;
				run;

				data _null_;
					set st;
					call symputx("lo&_I_" ,P_lo);
					call symputx("hi&_I_",P_hi);
				run;

				data vars2;
					set vars2;
					if name="&&iv&_I_" then do;
						P&p_lo=&&lo&_I_;
						P&p_hi=&&hi&_I_;
					end;
				run;
				%if %sysevalf(&&lo&_I_=&&hi&_I_) %then %put variable has
					constant value and will be excluded;
				%else %do;
					%pcont(&insdn,&&iv&_I_,&&miss&_I_);
				%end;
			%end;
		%end;

		%** Create permanent variable summary file **;
		data out.CE2_continuous_vars;
			set vars2;
		run;
		%** Write out list of variables to keep to code file **;
		%let rck=0;

		data _NULL_;
			set out.CE2_continuous_vars (where=(missing(new_var)=0)) end=eof;
			m=mod(_N_,7);
			length v $256.;
			retain v ;
			if m=1 then v=new_var;
			else v=strip(v) || " " || strip(new_var);
			FILE "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
			if _N_=1 then do;
				PUT ' ' '0d'x;
				PUT ' ' '0d'x;
				PUT '%LET KEEPLIST_C =';
			end;
			if m=0 or eof then PUT v '0d'x;
			if eof then PUT ';';
			if eof then call symputx("rck" ,_N_);
		run;
		%if &rck=0 %then %do;
			data _null_;
				FILE "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
				PUT ' ' '0d'x;
				PUT '%LET KEEPLIST_C = ;' '0d'x;
			run;
		%end;
		%put number of recoded continuous variables=&rck;

		%* Capture old and new variable name;
		data tmp (keep=orig_var orig_label variable);
			length variable orig_var $32.;
			length orig_label $256.;
			set out.CE2_continuous_vars (where=(missing(variable)=0)
				rename=(name=orig_var label=orig_label new_var=variable));
		run;

		data keep_vars;
			set keep_vars tmp;
		run;

		proc datasets nolist;
			delete vars vars2 tmp;
		quit;
		run;
	%end;
	%** If there are no continuous variables to process **;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Continuous_Var_Recode.txt" mod;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_C = ;' '0d'x;
		run;
	%end;
%mend cont_cntl;

** Master EDA and Recode macro **;
%macro CE_EDA_Recode(insdn=out.CE1_Resampled);
	%** Extract working portion of file **;
	data workfile;
		set &insdn (where=(mod_val_test ^= 3));
	run;

	%* Get global numbers;
	%global nobs overall_avg minbinn;

	proc sql noprint;
		select count(*) into : nobs from workfile;
		select mean(&dep_var) into : overall_avg from workfile;
	quit;

	data _null_;
		call symputx("minbinn",max(&minbinnc,(&nobs*&minbinnp)));
		** minimum bin size for nominal;
	run;

	%** Initialize profiling dataset **;
	%if %upcase(&profiling)=Y %then %do;
		data out.CE2_profile;
			set _NULL_;
			length type $10.;
			length variable $32.;
			length label $256.;
			length category $256.;
			format count comma8.;
			format percent percent8.2;
			%if %upcase(&binary_dv)=Y %then %do;
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
	data keep_vars;
		set _NULL_;
		length orig_var $32.;
		length orig_label $256.;
		length variable $32.;
	run;

	%if %symexist(binvar) %then %do;
		%bin_cntl(workfile);
	%end;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Binary_Var_Recode.txt" LRECL=256;
			PUT
				"*********************************************************************************;"
				'0d'x/
				"****                        CE2_BINARY_VAR_RECODE.TXT                        ****;"
				'0d'x/
				"*********************************************************************************;"
				'0d'x;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_B = ;' '0d'x;
		run;
	%end;

	%if %symexist(nomvar) %then %do;
		%nom_cntl(workfile);
	%end;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Nominal_Var_Recode.txt" LRECL=256;
			PUT
				"*********************************************************************************;"
				'0d'x/
				"****                        CE2_NOMINAL_VAR_RECODE.TXT                       ****;"
				'0d'x/
				"*********************************************************************************;"
				'0d'x;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_N = ;' '0d'x;
		run;
	%end;

	%if %symexist(ordvar) %then %do;
		%ord_cntl(workfile);
	%end;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Ordinal_Var_Recode.txt" LRECL=256;
			PUT
				"*********************************************************************************;"
				'0d'x/
				"****                       CE2_ORDINAL_VAR_RECODE.TXT                        ****;"
				'0d'x/
				"*********************************************************************************;"
				'0d'x;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_O = ;' '0d'x;
		run;
	%end;

	%if %symexist(contvar) %then %do;
		%cont_cntl(workfile);
	%end;
	%else %do;
		data _null_;
			FILE "&Path_output.CE2_Continuous_Var_Recode.txt" LRECL=256;
			PUT
				"*********************************************************************************;"
				'0d'x/
				"****                      CE2_CONTINUOUS_VAR_RECODE.TXT                      ****;"
				'0d'x/
				"*********************************************************************************;"
				'0d'x;
			PUT ' ' '0d'x;
			PUT '%LET KEEPLIST_C = ;' '0d'x;
		run;
	%end;

	%** Create dataset with recoded variables **;
	data out.CE2_Recoded;
		set &insdn;
		%inc "&Path_output.CE2_Binary_Var_Recode.txt";
		%inc "&Path_output.CE2_Nominal_Var_Recode.txt";
		%inc "&Path_output.CE2_Ordinal_Var_Recode.txt";
		%inc "&Path_output.CE2_Continuous_Var_Recode.txt";
		keep &keep_list &KEEPLIST_B &KEEPLIST_N &KEEPLIST_O &KEEPLIST_C;
	run;

	%** Finalize variable name lookup dataset **;
	proc contents data=out.CE2_Recoded out=vars (keep=name label) noprint;
	run;

	proc sql;
		create table out.CE2_vars as select a.variable, b.label, a.orig_var,
			a.orig_label from keep_vars a, vars b where a.variable=b.name order
			by a.variable;
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
		ods tagsets.excelxp options(sheet_name="Binary Variables"
			absolute_column_width="&l1,8,8,&l2,&l3,8,8" Frozen_Headers='Yes'
			Frozen_RowHeaders='1');

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
		ods tagsets.excelxp options(sheet_name="Nominal Variables"
			absolute_column_width="&l1,8,8,&l2,&l3,8,8,8" Frozen_Headers='Yes'
			Frozen_RowHeaders='1');

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
		ods tagsets.excelxp options(sheet_name="Ordinal Variables" %if
			%upcase(&stdmethodC) ^= NO %then %do;
		%if %upcase(&binary_dv)=Y %then %do;
		absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,8,8,9,8,8,8,8"
			%end;
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
		%if %upcase(&stdmethodO) ^= NO %then %do;
			var Location Scale / style={tagattr='format:0.00'};
		%end;
		var miss_impute / style={tagattr='format:0.0000'};
		%if %upcase(&Binary_dv)=Y %then %do;
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
		ods tagsets.excelxp options(sheet_name="Continuous Variables" %if
			%upcase(&stdmethodC) ^= NO %then %do;
		%if %upcase(&binary_dv)=Y %then %do;
		absolute_column_width="&l1,8,8,&l2,&l3,8,8,8,8,8,8,8,8,8,9,8,8,8,8"
			%end;
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
		var P&p_lo P&p_hi var_lb var_ub var_median /
			style={tagattr='format:0.00'};
		%if %upcase(&stdmethodC) ^= NO %then %do;
			var Location Scale / style={tagattr='format:0.00'};
		%end;
		var miss_impute / style={tagattr='format:0.0000'};
		%if %upcase(&Binary_dv)=Y %then %do;
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
		create table out.CE2_corr as select a.name as variable, a.label,
			b.&dep_var as correlation from vars a, corr b where a.name=b._name_
			order by abs(b.&dep_var) desc;
	quit;

	proc sql noprint;
		select max(length(variable)) into: l1 from out.CE2_corr;
		select max(length(label)) into: l2 from out.CE2_corr;
	quit;

	ods listing close;
	ods Tagsets.ExcelxP body="&Path_output.CE2_Correlation_report.xls"
		style=sasweb;
	ods tagsets.excelxp options(sheet_name="Correlations"
		absolute_column_width="&l1,&l2,8" Frozen_Headers='Yes'
		Frozen_RowHeaders='1');

	proc print data=out.CE2_corr noobs;
		var variable label;
		var correlation / style={tagattr='format:0.0000'};
	run;
	ods Tagsets.ExcelxP close;
	ods listing;

	%if %upcase(&profiling)=Y %then %do;
		proc sql noprint;
			select max(length(variable)) into: l1 from out.CE2_profile;
			select max(length(label)) into: l2 from out.CE2_profile;
			select max(length(category)) into: l3 from out.CE2_profile;
		quit;

		%** Create Profile report in Excel;

		%* Order file;
		proc sql;
			create table ord as select variable, max(abs(index))-min(abs(index))
				as diff from out.CE2_profile group by variable order by diff
				desc;
		quit;

		data ord;
			set ord;
			sord=_n_;
		run;

		data prof;
			set out.CE2_profile;
			s2=_n_;
		run;

		proc sql;
			create table prof2 as select a.* from prof a, ord b where a.variable
				=b.variable order by b.sord, a.s2;
		quit;

		data out.CE2_profile;
			set prof2 (drop=s2);
		run;

		proc datasets nolist;
			delete prof prof2;
		quit;
		run;

		ods listing close;
		ods Tagsets.ExcelxP body="&Path_output.CE2_Profile_report.xls"
			style=sasweb;

		ods tagsets.excelxp options(sheet_name="Binary Variables"
			absolute_column_width="&l1,&l2,&l3,8,8,9,8,6" Frozen_Headers='Yes'
			Frozen_RowHeaders='1');

		proc print data=out.CE2_profile noobs;
			where type="Binary";
			var variable label category;
			var count / style={tagattr='format:#,##0'};
			var percent Average_DV index star;
		run;
		ods tagsets.excelxp options(sheet_name="Nominal Variables"
			absolute_column_width="&l1,&l2,&l3,8,8,9,8,6" Frozen_Headers='Yes'
			Frozen_RowHeaders='1');

		proc print data=out.CE2_profile noobs;
			where type="Nominal";
			var variable label category;
			var count / style={tagattr='format:#,##0'};
			var percent Average_DV index star;
		run;
		ods tagsets.excelxp options(sheet_name="Ordinal Variables"
			absolute_column_width="&l1,&l2,&l3,8,8,9,8,6" Frozen_Headers='Yes'
			Frozen_RowHeaders='1');

		proc print data=out.CE2_profile noobs;
			where type="Ordinal";
			var variable label category;
			var count / style={tagattr='format:#,##0'};
			var percent Average_DV index star;
		run;
		ods tagsets.excelxp options(sheet_name="Continuous Variables"
			absolute_column_width="&l1,&l2,&l3,8,8,9,8,6" Frozen_Headers='Yes'
			Frozen_RowHeaders='1');

		proc print data=out.CE2_profile noobs;
			where type="Continuous";
			var variable label category;
			var count / style={tagattr='format:#,##0'};
			var percent Average_DV index star;
		run;

		ods Tagsets.ExcelxP close;
		ods listing;
	%end;

	proc datasets nolist;
		delete workfile corr vars keep_vars;
	quit;
	run;
%MEND;
