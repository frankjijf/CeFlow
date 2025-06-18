
import pandas as pd
import numpy as np
import os

# Loading variable lists
import Titanic_test.contvar_list as contvar_list
import Titanic_test.nomvar_list as nomvar_list
import Titanic_test.binvar_list as binvar_list
import Titanic_test.ordvar_list as ordvar_list

# Assigning paths and dataset names
lib_dir = r"D:\OneDrive\CE_PROJECT\Python_CE_Project_GH\SAS2PYTHON\CE_PY_CODE\Titanic_test"
out_dir = os.path.join(lib_dir, "ce_output")
os.makedirs(lib_dir, exist_ok=True)

# Create output directory if it doesn't exist
log_dir = os.path.join(lib_dir, "log")
os.makedirs(log_dir, exist_ok=True)
log_path_01 = os.path.join(log_dir, "01_CE_Sampling_Log_File.log")
lst_path_01 = os.path.join(log_dir, "01_CE_Sampling_LST_File.lst")
log_path_02 = os.path.join(log_dir, "02_CE_EDA_Recode_Log_File.log")
lst_path_02 = os.path.join(log_dir, "02_CE_EDA_Recode_LST_File.lst")

# Load the dataset
inds = pd.read_csv(os.path.join(lib_dir, "Titanic-Dataset.csv"))  # Your SAS input dataset name (assuming CSV)
inds['IsFemale'] = (inds['Sex']=='female').astype(int)
inds['IsChild'] = (inds['Age'] < 18).astype(int)

# Mandatory Macro Variables: These must be set
config = {
    # Variable Lists
        "contvar" : contvar_list.contvar_list(),  # Your continuous variable list
        "nomvar" : nomvar_list.nomvar_list(),  # Your nominal variable list
        "binvar" : binvar_list.binvar_list(),  # Your binary variable list
        "ordvar" : ordvar_list.ordvar_list(),  # Your ordinal variable list
    # General Macro Variables
        "path_output" : out_dir,  # your output destination folder
        "inds" : inds,  # your input dataset name
        "id" : 'PassengerId',  # Unique ID
        "dep_var" : 'Survived',  # Your dependent variable
        "binary_dv" : 'Y',  # Your dependent variable type (Y/N)
        "weight" : None,  # If you want to use weights, set this variable to the name of your weight variable
    # Macro 1: Sampling macro variables
        "split_portion": 0.7,  # Your split portion for modeling portion
        "exclusion_if": "inds[dep_var].isna()",  # Data exclusion condition; set an actual condition if needed, e.g., "data[dep_var].notna()"
        "DS_present": 'N',  # Change to 'Y' only if your data contains the standard DS modeling bundle
    # Optional Macro Variables: These all have defaults that can be used
    # General Macro Variables
        "prefix" : 'R1_',  # Your variable recoding prefix
        "keep_list" : None,  # Add any additional variables you want to keep for analysis purposes
    # Macro 1: Sampling macro variables
        "path_DS": "/mnt/projects/shared/pst_qmgisi/Modeling/CE/",  # Location of standard DS recodes. Do not change
        "bootstrap": 'Y',  # Request oversampling, bootstrap if responders size is small (Y/N). Binary DV only
        "oversampled_rr": 0.1,  # Your desired oversample response rate on the modeling part of data, 0.1-0.5
        "min_num_resp": 2000,  # Minimum number of responders if sampling
        "seed": 123456,  # Random selection seed if sampling
    # Macro 2: Recoding macro variables
        "profiling": 'Y',  # Request profiling report on all variables (Y/N)
        "missrate": 0.75,  # Maximum missing rate allowed
        "concrate": 0.9,  # Maximum amount of file that can be in a single value of an independent variable
        "valcnt": 50,  # Maximum number of unique values allowed. Set to 0 to allow any number
        "minbinnc": 500,  # Minimum count in a bin to be usable. Set to 0 to use minbinpct only
        "minbinnp": 0.05,  # Minimum percent of file in a bin to be usable. Set to 0 to use minbincnt only
        "talpha": 0.05,  # T-Test significance level for collapse of bins
        "bonfer": 'N',  # Do Bonferroni adjustment for talpha? (Y/N)
        "nom_method": 'INDEX',  # Recoding method for nominal variables: Binary, Index or Mean
        "pvalue": 0.05,  # P-value threshold to include variable
        "min_size": 500,  # Minimum missing group size for Equal Response imputation
        "num_category": 10,  # Maximum number of categories for profiling variables
        "equal_dist": 'N',  # Use equal distance for dividing variables into groups for profiling (Y/N)
        "p_lo": 1,  # Lower percentile for checking constant value
        "p_hi": 99,  # Upper percentile for checking constant value
        "impmethodC": 'median',  # What method to use for missing imputation for continuous variables?
        # Options are ER for Equal Response or any proc stdize method
        "stdmethodC": 'STD',  # Standardization options: any method allowed in proc stdize or NO to skip
        "cap_flrC": 'Y',  # Do you want to do capping/flooring to handle outliers? (Y/N)
        "transformationC": 'Y',  # Include transformed variables in evaluation (Y/N)
        "impmethodO": 'mean',  # What method to use for missing imputation for ordinal variables?
        "stdmethodO": 'No',  # Standardization options: any method allowed in proc stdize or NO to skip
        "cap_flrO": 'N',  # Do you want to do capping/flooring to handle outliers? (Y/N)
        "transformationO": 'N',  # Include transformed variables in evaluation? (Y/N)
}

################# RUN #################

from CE1_sampling import CE_Sampling
from CE_log import printto
from CE2_1_pnum_prof123 import pnum, prof1, prof2, prof3
from CE2_2_pbin_bin_cntl import bin_cntl
from CE2_3_pnorm_nom_cntl import nom_cntl
from CE2_4_pord_ord_cntl import ord_cntl
from CE2_5_pcon_cont_cntl import cont_cntl
from CE2_6_CE_EDA_Recode import CE_EDA_Recode

# 调用 CE_Sampling 函数时传递这些参数
with printto(log=log_path_01, lst=lst_path_01) as logger:
    CE1_Resampled, CE1_Sample_Rate = CE_Sampling(
        config  = config,
        logger  = logger)

with printto(log=log_path_02, lst=log_path_02) as logger:
    CE2_Recoded, profile_df, var_lookup_df = CE_EDA_Recode(
        indsn   = CE1_Resampled,
        config  = config,
        logger  = logger
    )


















# Macro 3: Variable Reduction macro variables
fast_opt = 'N'  # Fast option will turn off all tests except multivariate regression
# Optional Macro Variables: These all have defaults that can be used
samplesize = 50000  # Sample size to use for variable reduction
redu_weight = 'N'  # Use weights in variable reduction? (Y/N)
sources = 3  # Minimum number of sources to be selected
# Univariate regression option
univ_reg = 'Y'  # Use univariate regression to choose variables? (Y/N)
maxpuni = 0.0001  # Maximum p-value correlation for selecting via univariate regression
# Correlation option
correlation = 'Y'  # Use correlation to choose variables? (Y/N)
corrcut = 0.01  # Minimum correlation between independent variable and dependent variable
# Principal components option
principal = 'Y'  # Use principal components to choose variables? (Y/N)
nprin = 10  # Number of principal components desired
minprin = 0.5  # Minimum factor correlation for selecting via principal component
# Cluster option
cluster = 'Y'  # Use cluster analysis to choose variables? (Y/N)
maxc = 20  # Number of clusters desired
maxratio = 0.5  # Maximum R-squared ratio for selecting variables via clustering
# Linear regression option
regression = 'Y'  # Use linear regression to choose variables? (Y/N)
alphareg = 0.05  # Alpha level for forward selection in linear regression
# Logistic regression option - only applicable if binary dependent variable
logistic = 'Y'  # Use logistic regression to choose variables? (Y/N)
alphalog = 0.05  # Alpha level for forward selection in logistic regression
# Information value option
information = 'Y'  # Use information value to choose variables? (Y/N)
decile = 20  # Number of groups to use when calculating information values
infvcut = 0.01  # Minimum information value for selecting via information value
# Maximum correlation between independent variables option
ind_correlation = 'Y'  # Exclude variables with high correlation to others? (Y/N)
maxcorr = 0.7  # Maximum correlation allowed between independent variables
# Maximum correlation to dependent variable option
ind_dv_corr = 'Y'  # Exclude variables with high correlation to dependent variable? (Y/N)
max_dv_corr = 0.7  # Maximum correlation allowed to dependent variable

# Macro 4: Model selection and tuning
# Optional Macro Variables: These all have defaults that can be used
sel_alpha = 0.05  # Alpha level when selecting/removing variables into the model
includelist = []  # List of variables that have to be in the final model
startlist = []  # List of variables that need to be in the first step
excludelist = []  # List of variables to be excluded from modeling
criteria = 'c'  # Metric to use during variable tuning. Default is c for binary, AdjRsq for other
threshold = 0  # Minimum change in evaluation metric to include variable
SQL_join = 'union'  # Type of join to use between file portions during variable tuning
minimp = 0.01  # Minimum relative importance to keep variable
graph_plot = 'Y'  # Include graphing? (Y/N)

# Macro 5: Final Model selection
# Optional Macro Variables: These all have defaults that can be used
fin_alpha = 0.05  # Alpha level when selecting/removing variables into the model
method = 'stepwise'  # Model build method
fin_num_category = 10  # Maximum number of categories for profiling variables
fin_equal_dist = 'N'  # Use equal distance for dividing variables into groups for profiling (Y/N)









# 2. EDA, profiling and Recode
LogFile = path_output + "02_CE_Var_EDA_Recode_Log_File.log"
LstFile = path_output + "02_CE_Var_EDA_Recode_List_File.lst"
# Placeholder for EDA and recode logic
CE2_Recoded = CE1_Resampled.copy()  # Example, modify as needed

# 3. Variable reduction and ranking
LogFile = path_output + "03_Var_Redu_Log_File.log"
LstFile = path_output + "03_Var_Redu_List_File.lst"
# Placeholder for variable reduction logic
varlist_redu = ['placeholder_var1', 'placeholder_var2']  # Example, replace with actual logic
# Implement actual variable reduction logic here

# 4. Model selection and tuning
LogFile = path_output + "04_Model_Selection_Log_File.log"
LstFile = path_output + "04_Model_Selection_List_File.lst"
# Model selection logic (example: logistic regression for binary DV)
model = LogisticRegression()
X = CE2_Recoded[varlist_redu]
y = CE2_Recoded[dep_var]
model.fit(X, y)

# 5. Final Model build and validation on test sample
LogFile = path_output + "05_Final_Model_Fit_Validation_Log_File.log"
LstFile = path_output + "05_Final_Model_Fit_Validation_LST_File.lst"
# Placeholder for final model validation logic
varlist_final = ['placeholder_final_var1', 'placeholder_final_var2']  # Example, replace with actual logic
X_final = CE2_Recoded[varlist_final]
y_final = CE2_Recoded[dep_var]
final_model = LogisticRegression()
final_model.fit(X_final, y_final)

# Example output
output_predictions = final_model.predict(X_final)

