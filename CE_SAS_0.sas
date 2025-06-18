
%let loc=/mnt/projects/public/offshore_office/frankji/CE exercise/exercise/sasdata/; ** Project folder;

libname lib   "&loc";                                      ** SAS dataset libname;
libname out   "&loc.ce_output";                                    ** CE output libname;

* Comment out any variable types that are not present in the data;
%include "&loc.contvar_list.txt";                                  ** List of continuous variables to include in analysis;
%include "&loc.ordvar_list.txt";                                   ** List of ordinal variables to include in analysis;
%include "&loc.nomvar_list.txt";                                   ** List of nominal variables to include in analysis;
%include "&loc.binvar_list.txt";                                   ** List of binary variables to include in analysis;

***  The following are macro variables needed for CE  ***;

*** Mandatory Macro Variables: These must be set ***;
* General Macro Variables *;
%let path_output=&loc.ce_output/;    ** your output destination folder;
%let inds=lib.randomsample;               ** your sas input dataset name;
%let id = rid;                      ** Unique ID;
%let dep_var=dv_dpa_purchase_flag_in_attr;                   ** your dependant variable;
%let binary_dv=Y;                     ** your dependent variable type (Y/N);
%let weight= ;                        ** If you want to use weights, set this variable to the name of your weight variable;

*** Macro 1: Sampling macro variables ***;
%let split_if = if ranuni(12345)<0.7; ** your split condition for modeling portion;
%let exclusion_if=;     ** Data exclusion if condition applies. e.g. IF &Dep_var = .;
%let DS_Present= N;                   ** change to Y only if your data contain standard DS modeling bundle;

*** Macro 3: Variable Reduction macro variables ***;
%let fast_opt = N;                    ** Fast option will turn off all tests except multivariate regression

*** Optional Macro Variables: These all have defaults that can be used ***;
* General Macro Variables *;
%let prefix=R1_;                      ** your variable recoding prefix;
%let keep_list=&id &dep_var &weight mod_val_test;  ** Add any additional variables you want to keep for analysis purposes;









*** Macro 1: Sampling macro variables ***;
%let path_DS=/mnt/projects/shared/pst_qmgisi/Modeling/CE/;  ** Location of standard DS recodes. Do Not Change;
%let bootstrap= N;                    ** request oversamping, bootstrap if responders size is small (Y/N). Binary dv only;
%let oversampled_rr=0.05;             ** your desired oversample response rate on the modeling part of data, 0.1-0.5;
%let min_num_resp=500;                ** Minimum number of responders if sampling;
%let seed=123456;                     ** Random selection seed if sampling;

*** Macro 2: Recoding macro variables ***;
%let Profiling = Y;                   ** Request profiling report on all variables (Y/N);
%let missrate = .75;                  ** maximum missing rate allowed;
* Binary and Nominal variable types;
%let concrate = .9;                   ** Maximum amount of file that can be in a single value of an independent variable;
* Nominal variable types;
%let valcnt = 50;                     ** Maximum number of unique values allowed. Set to 0 to allow any number;
%let minbinnc = 500;                  ** Minimum count in a bin to be usable. set to 0 to use minbinpct only;
%let minbinnp = .05;                  ** Minimum percent of file in a bin to be usable. set to 0 to use minbincnt only;
%let talpha  = 0.05;                  ** T Test significance level for collapse of bins;
%let bonfer = N;                      ** Do Bonferoni adjustment for talpha? (Y/N);
%let nom_method = INDEX;              ** Recoding method for nominal variables: Binary, Index or Mean;
* Continuous and Ordinal variable types;
%let pvalue = .05;                    ** Pvalue threshold to include variable;
%let min_size = 500;                  ** Minimum missing group size for Equal Response imputation;
%let num_category = 10;               ** Maximum number of categories for profiling variables;
%let equal_dist = N;                  ** Use equal distance for dividing variables into groups for profiling (Y/N);
* Continuous variable types;
%let p_lo =  1;                       ** Lower percentile for checking constant value;
%let p_hi = 99;                       ** Upper percentile for checking constant value;
%let impmethodC = median;               ** What method to use for missing imputation for continuous variables?;
                                      ** Options are ER for Equal Reponse or any proc stdize method;
         ** The location associated with the method will be used for missing values. Be sure to include c or p value if needed;
         **    Method      Location 
         **    -------     ---------
         **    MEAN        Mean
         **    MEDIAN      Median
         **    SUM         0
         **    EUCLEN      0
         **    USTD        0
         **    STD         Mean
         **    RANGE       Minimum
         **    MIDRANGE    Midrange
         **    MAXABS      0
         **    IQR         Median
         **    MAD         Median
         **    ABW(c)      Biweight one-step M-estimate
         **    AHUBER(c)   Huber one-step M-estimate
         **    AWAVE(c)    Wave one-step M-estimate
         **    AGK(p)      Mean
         **    SPACING(p)  Mid-minimum spacing
         **    L(p)        L(p)  ;
%let stdmethodC = STD;                ** Standardization options: any method allowed in proc stdize or NO to skip;
%let cap_flrC = Y;                    ** Do you want to do capping / flooring to handle outliers? (Y/N);
%let transformationC = Y;             ** Include transformed variables in evaluation (Y/N);
* Ordinal variable types;
%let impmethodO = mean;               ** What method to use for missing imputation for continuous variables?;
                                      ** Options are ER for Equal Reponse or any proc stdize method;
%let stdmethodO = No;                 ** Standardization options: any method allowed in proc stdize or NO to skip;
%let cap_flrO = N;                    ** Do you want to do capping / flooring to handle outliers? (Y/N);
%let transformationO = N;             ** Include transformed variables in evaluation? (Y/N);

*** Macro 3: Variable Reduction macro variables ***;
%let samplesize = 50000;              ** Sample size to use for variable reduction;
%let redu_weight = N;                 ** Use weights in variables reduction? (Y/N);
%let sources = 3;                     ** Minimum number of sources to be selected;
* Univariate regression option;
%let univ_reg = Y;                    ** Use univariate regression to choose variables? (Y/N);
%let maxpuni = .0001;                 ** Maximum p value correlation for selecting via univariate regression;
* Correlation option;
%let correlation = Y;                 ** Use correlation to choose variables? (Y/N);
%let corrcut = .01;                   ** Minimum correlation between independent variable and dependent variable;
* Principal components option;
%let principal = Y;                   ** Use principal components to choose variables? (Y/N);
%let nprin = 10;                      ** Number of principal components desired;
%let minprin = .5;                    ** Minimum factor correlation for selecting via principal component;
* Cluster option;
%let cluster = Y;                     ** Use cluster analysis to choose variables? (Y/N);
%let maxc = 20;                       ** Number of clusters desired;
%let maxratio = .5;                   ** Maximum R Squared ratio for selecting variables via clustering;
* Linear regression option;
%let regression = Y;                  ** Use linear regression to choose variables? (Y/N);
%let alphareg = .05;                  ** Alpha level for forward selection in linear regression;
* Logistic regression option - only applicable if binary dependent variable;
%let logistic = Y;                    ** Use logistic regression to choose variables? (Y/N);
%let alphalog = .05;                  ** Alpha level for forward selection in logistic regression;
* Information value option;
%let information = Y;                 ** Use information value to choose variables? (Y/N);
%let decile = 20;                     ** Number of groups to use when calculate information values;
%let infvcut = .01;                   ** Minimum information value for selecting via information value;
* Maximum correlation between independent variables option;
%let ind_correlation = Y;             ** Exclude variables with high correlation to others? (Y/N);
%let maxcorr = .7;                    ** Maximum correlation allowed between independent variables;
* Maximum correlation to dependent variable option;
%let ind_dv_corr=Y;                   ** Exclude variables with high correlation to dependent variable? (Y/N);
%let max_dv_corr=.7;                  ** Maximum correlation allowed to dependent variable;

*** Macro 4: Model selection and tuning ***;
%let sel_alpha = .05;                 ** Alpha level when selecting/removing variable into the model;
%let includelist = ;                  ** List of variables that has to be in the final model;
%let startlist = ;                    ** List of variables that need to be in the first step;
%let excludelist = ;                  ** List of variables to be excluded from modeling;
%let criteria = c;                    ** Metric to use during variable tuning. Default is c for binary, AdjRsq for other;
%let threshold = 0;                   ** Minimum change in evaluation metric to include variable;
%let SQL_join = union;                ** Type of join to use between file portions during variable tuning;
%let minimp = .01;                    ** Minimum relative importance to keep variable;
%let graph_plot = Y;                  ** Include graphing? (Y/N);

*** Macro 5: Final Model selection ***;
%let fin_alpha = .05;                 ** Alpha level when selecting/removing variable into the model;
%let method = stepwise;               ** Model build method;
%let fin_num_category = 10;           ** Maximum number of categories for profiling variables;
%let fin_equal_dist = N;              ** Use equal distance for dividing variables into groups for profiling (Y/N);

ods html close;
ods listing;

%inc "/mnt/projects/shared/pst_qmgisi/Modeling/CE/CE_Macros_v5.sas";


*** 1.Sampling ***;
%let LogFile = "&path_output.01_CE_Sampling_Log_File.log";
%let LstFile = "&path_output.01_CE_Sampling_LST_File.lst";
proc printto log=&LogFile new; run;
proc printto print=&LstFile new; run;

%CE_Sampling(inds=&inds, outds=out.CE1_Resampled)

proc printto;
run;

*** 2.EDA, profiling and Recode ***;
%let LogFile = "&path_output.02_CE_Var_EDA_Recode_Log_File.log";
%let LstFile = "&path_output.02_CE_Var_EDA_Recode_List_File.lst";
proc printto log=&LogFile new; run;
proc printto print=&LstFile new; run;

%CE_EDA_RECODE(INSDN=out.CE1_Resampled);

proc printto;
run;

*** 3.Variable reduction and ranking ***;
%let LogFile = "&path_output.03_Var_Redu_Log_File.log";
%let LstFile = "&path_output.03_Var_Redu_List_File.lst";
proc printto log=&LogFile new; run;
proc printto print=&LstFile new; run;

%CE_Var_Redu(insdn=out.CE2_Recoded);

proc printto;
run;

*** 4.Model selection and tuning ***;
%let LogFile = "&path_output.04_Model_Selection_Log_File.log";
%let LstFile = "&path_output.04_Model_Selection_List_File.lst";
proc printto log=&LogFile new; run;
proc printto print=&LstFile new; run;

%inc "&path_output.CE3_Varlist_redu.txt";
%CE_Model_Val(out.CE2_Recoded, &varlist_redu);

proc printto;
run;

*** 5.Final Model build and validation on test sample ***;
%let LogFile = "&path_output.05_Final_Model_Fit_Validation_Log_File.log";
%let LstFile = "&path_output.05_Final_Model_Fit_Validation_LST_File.lst";
proc printto log=&LogFile new; run;
proc printto print=&LstFile new; run;

%inc "&path_output.CE4_Varlist_Final.txt";
%CE_Model_Lift(insdn=out.CE2_Recoded, varlist=&varlist_final)

proc printto;
run;
