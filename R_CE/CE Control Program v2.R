## ***********************************************************************;
## ***********************************************************************;
## ****                 CONTENT EVALUATOR - R                         ****;
## ****   VERSION 2.0:  UPDATED ON 9/12/2014 Nan Flaaten              ****;
## ***********************************************************************;
## ***********************************************************************;

# Load required packages
require(XLConnect)
require(plyr)
require(subselect)
require(car)

# Load CE functions
source("/mnt/projects/shared/pst_qmgisi/Modeling/RCE/CE Functions v2.R", echo = F)

# Set working directory to project directory (where output will be written out to)
# Make sure your data lists are in this location
setwd("/mnt/projects/locked/pst_qmgisi/Modeling/CE Rebuild/binary/R_work")

## Step 1: Sampling
lib <- "/mnt/projects/locked/pst_qmgisi/Modeling/CE Rebuild/sasdata" # Path to input file
f.sample(dat = file.path(lib, "bin_file.xdf"),        # your input dataset name - xdf
         outds = "CE1_Resampled.xdf",                 # output dataset name
         dep_var = "resp",                            # your dependant variable
         binary_dv = "Y",                             # your dependent variable type (Y/N)
         include_if = expression(!is.na(resp)),       # data inclusion if condition applies - either expression(..) or NULL
         split_if = .7,                               # your split condition for modeling portion - either expression are decimal
         ds_present = "Y",                            # change to Y if your data contains standard DS modeling bundle
         bootstrap = "N",                             # request oversamping, bootstrap if responders size is small (Y/N). Binary dv only
         oversampled_rr = .20,                        # your desired oversample response rate on the modeling part of data, 0.1-0.5
         min_num_resp = 1000,                         # minimum number of responders if sampling
         seed = 123456)                               # random selection seed if sampling

## Step 2: EDA and variable recoding
f.recode(dat = "CE1_Resampled.xdf",                   # your input dataset name - output from step 1
         keep_list=c("KL_id","resp", "wgt", "mod_val_test"),   # variables you want to keep for analysis purposes - 
                                                      # make sure you include your id variable, dependent variable and mod_val_test
         dep_var = "resp",                            # your dependant variable
         binary_dv = "Y",                             # your dependent variable type (Y/N)
         prefix = "R1_",                              # your variable recoding prefix
         missrate = .75,                              # maximum missing rate allowed
         concrate = .9,                               # maximum amount of file that can be in a single value of an independent variable
         valcnt = 50,                                 # maximum number of unique values allowed. Set to 0 to allow any number
         minbinnc = 500,                              # minimum count in a bin to be usable. set to 0 to use minbinpct only
         minbinnp = .05,                              # Minimum percent of file in a bin to be usable. set to 0 to use minbincnt only
         talpha  = .05,                               # T Test significance level for collapse of bins
         nom_method = "INDEX",                        # recoding method for nominal variables: Binary, Index or Mean
         min_size = 500,                              # minimum missing group size for Equal Response imputation
         cap_flrO = "Y",                              # ordinal - do capping / flooring to handle outliers? (Y/N)
         transformationO = "N",                       # ordinal - include transformed variables in evaluation? (Y/N)
         impmethodO = "MEDIAN",                       # ordinal - what method to use for missing imputation?
         stdmethodO = "NO",                           # ordinal - standardization method? (NO to skip)
         p_lo = .01,                                  # continuous - lower percentile for checking constant value
         p_hi = .99,                                  # continuous - upper percentile for checking constant value
         cap_flrC = "Y",                              # continuous - do capping / flooring to handle outliers? (Y/N)
         transformationC = "Y",                       # continuous - include transformed variables in evaluation? (Y/N)
         impmethodC = "MEAN",                         # continuous - what method to use for missing imputation?
         stdmethodC = "STD",                          # continuous - standardization method? (NO to skip)
         profiling = "Y",                             # request profiling report on all variables (Y/N)
         equal_dist = "N",                            # use equal distance for dividing variables into groups for profiling (Y/N)
         num_category = 10)                           # maximum number of categories for profiling variables

## Step 3: Variable reduction
f.var_redu(dat = "CE2_Recoded.xdf",                   # your input dataset name - output from step 2
           dep_var = "resp",                          # your dependant variable
           binary_dv = "Y",                           # your dependent variable type (Y/N)
           samplesize = 50000,                        # sample size to use for variable reduction
           redu_weight = "N",                         # use weights in variables reduction? (Y/N)
           weight = "wgt",                            # name of weight variable
           sources = 3,                               # minimum number of sources to be selected
           maxnum = 1000,                             # maximum number of variables to keep
           maxcorr = .7,                              # maximum correlation allowed between independent variables
           ind_dv_corr = "Y",                         # exclude variables with high correlation to dependent variable? (Y/N)
           max_dv_corr = .7,                          # maximum correlation allowed to dependent variable
           univ_reg = "N",                            # use univariate regression to choose variables? (Y/N)
           maxpuni = .05,                             # maximum p value correlation for selecting via univariate regression
           correlation = "Y",                         # use correlation to choose variables? (Y/N)
           corrcut = .01,                             # minimum correlation between independent variable and dependent variable
           factor = "Y",                              # use factor analysis to choose variables? (Y/N)
           nfact = 20,                                # number of factors desired
           minfact = .5,                              # minimum factor loading for selection
           regression = "Y",                          # use linear regression to choose variables? (Y/N)
           alphareg = .05,                            # alpha level for forward selection in linear regression
           logistic = "Y",                            # use logistic regression to choose variables? (Y/N) [binary dv only]
           alphalog = .05,                            # alpha level for forward selection in logistic regression [binary dv only]
           information = "Y",                         # use information value to choose variables? (Y/N)
           decile = 20,                               # number of groups to use when calculate information values
           infvcut = .01)                             # minimum information value for selecting via information value

## Step 4: Model selection and tuning
source("CE3_Varlist_redu.txt")  # Load variable list generated in step 3
f.model_val(insdn = "CE2_Recoded.xdf",                # your input dataset name - output from step 2
            varlist = varlist_redu,                   # list of variables to consider
            dep_var = "resp",                         # your dependant variable
            binary_dv = "Y",                          # your dependent variable type (Y/N)
            weight = NULL,                            # name of weight variable or NULL to not use weights
            sel_alpha = .05,                          # alpha level when selecting/removing variable into the model
            refit = F,                                # refit each step (T/F) [binary dv only]
            includelist = NA,                         # list of variables that has to be in the final model
            startlist = NA,                           # list of variables that need to be in the first step
            excludelist = NA,                         # list of variables to be excluded from modeling
            criteria = "c",                           # metric to use during variable tuning.
                                                      # Default is c for binary, AdjRsq for other
            threshold = 0,                            # minimum change in evaluation metric to include variable
            SQL_join = "union",                       # type of join to use between file portions during variable tuning
            minimp = .01,                             # minimum relative importance to keep variable
            graph_plot = "Y")                         # include graphing? (Y/N)

## Step 5: Final Model build and validation on test sample
source("CE4_Varlist_Final.txt")  # Load variable list generated in step 4
f.model_lift(insdn = "CE2_Recoded.xdf",               # your input dataset name - output from step 2
             varlist = varlist_final,                 # list of variables to consider
             dep_var = "resp",                        # your dependant variable
             binary_dv = "Y",                         # your dependent variable type (Y/N)
             weight = NULL,                           # name of weight variable or NULL to not use weights
             fin_alpha = .05,                         # alpha level when selecting/removing variable into the model
             refit = F,                               # refit each step (T/F) [binary dv only]
             method = "stepwise",                     # model build method
             fin_num_category = 10,                   # maximum number of categories for profiling variables
             fin_equal_dist = "N")                    # use equal distance for dividing variables into groups for profiling (Y/N)
                   