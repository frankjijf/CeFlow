import pandas as pd
import numpy as np
import os

# Loading variable lists
import Titanic_test.contvar_list as contvar_list
import Titanic_test.nomvar_list as nomvar_list
import Titanic_test.binvar_list as binvar_list
import Titanic_test.ordvar_list as ordvar_list

from ce1_sampling import CE_Sampling
from ce_log_tool import printto
from ce2_eda_recode import (
    pnum,
    prof1,
    prof2,
    prof3,
    bin_cntl,
    nom_cntl,
    ord_cntl,
    cont_cntl,
    CE_EDA_Recode,
)

# Assigning paths and dataset names
lib_dir = (
    r"D:\OneDrive\CE_PROJECT\Python_CE_Project_GH\SAS2PYTHON\CE_PY_CODE\Titanic_test"
)
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
inds = pd.read_csv(
    os.path.join(lib_dir, "Titanic-Dataset.csv")
)  # Your input dataset name (assuming CSV)
inds["IsFemale"] = (inds["Sex"] == "female").astype(int)
inds["IsChild"] = (inds["Age"] < 18).astype(int)

# Mandatory Macro Variables: These must be set
config = {
    # Variable Lists
    "cont_vars": contvar_list.contvar_list(),  # Your continuous variable list
    "nom_vars": nomvar_list.nomvar_list(),  # Your nominal variable list
    "bin_vars": binvar_list.binvar_list(),  # Your binary variable list
    "ord_vars": ordvar_list.ordvar_list(),  # Your ordinal variable list
    # General Macro Variables
    "path_output": out_dir,  # your output destination folder
    "inds": inds,  # your input dataset name
    "id": "PassengerId",  # Unique ID
    "dep_var": "Survived",  # Your dependent variable
    "binary_dv": "Y",  # Your dependent variable type (Y/N)
    "weight": None,  # If you want to use weights, set this variable to the name of your weight variable
    # Macro 1: Sampling macro variables
    "split_portion": 0.7,  # Your split portion for modeling portion
    "exclusion_if": "inds[dep_var].isna()",  # Data exclusion condition; set an actual condition if needed, e.g., "inds[dep_var].isna()"
    "DS_present": "N",  # Change to 'Y' only if your data contains the standard DS modeling bundle
    # Optional Macro Variables: These all have defaults that can be used
    # General Macro Variables
    "prefix": "R1_",  # Your variable recoding prefix
    "keep_list": None,  # Add any additional variables you want to keep for analysis purposes
    # Macro 1: Sampling macro variables
    "path_DS": "/mnt/projects/shared/pst_qmgisi/Modeling/CE/",  # Location of standard DS recodes. Do not change
    "bootstrap": "Y",  # Request oversampling, bootstrap if responders size is small (Y/N). Binary DV only
    "oversampled_rr": 0.1,  # Your desired oversample response rate on the modeling part of data, 0.1-0.5
    "min_num_resp": 2000,  # Minimum number of responders if sampling
    "seed": 123456,  # Random selection seed if sampling
    # Macro 2: Recoding macro variables
    "profiling": "Y",  # Request profiling report on all variables (Y/N)
    "missrate": 0.75,  # Maximum missing rate allowed
    "concrate": 0.9,  # Maximum amount of file that can be in a single value of an independent variable
    "valcnt": 50,  # Maximum number of unique values allowed. Set to 0 to allow any number
    "minbinnc": 500,  # Minimum count in a bin to be usable. Set to 0 to use minbinpct only
    "minbinnp": 0.05,  # Minimum percent of file in a bin to be usable. Set to 0 to use minbincnt only
    "talpha": 0.05,  # T-Test significance level for collapse of bins
    "bonfer": "N",  # Do Bonferroni adjustment for talpha? (Y/N)
    "nom_method": "INDEX",  # Recoding method for nominal variables: Binary, Index or Mean
    "pvalue": 0.05,  # P-value threshold to include variable
    "min_size": 500,  # Minimum missing group size for Equal Response imputation
    "num_category": 10,  # Maximum number of categories for profiling variables
    "equal_dist": "N",  # Use equal distance for dividing variables into groups for profiling (Y/N)
    "p_lo": 1,  # Lower percentile for checking constant value
    "p_hi": 99,  # Upper percentile for checking constant value
    "impmethodC": "median",  # What method to use for missing imputation for continuous variables?
    # Options are ER for Equal Response or any proc stdize method
    "stdmethodC": "STD",  # Standardization options: any method allowed in proc stdize or NO to skip
    "cap_flrC": "Y",  # Do you want to do capping/flooring to handle outliers? (Y/N)
    "transformationC": "Y",  # Include transformed variables in evaluation (Y/N)
    "impmethodO": "mean",  # What method to use for missing imputation for ordinal variables?
    "stdmethodO": "No",  # Standardization options: any method allowed in proc stdize or NO to skip
    "cap_flrO": "N",  # Do you want to do capping/flooring to handle outliers? (Y/N)
    "transformationO": "N",  # Include transformed variables in evaluation? (Y/N)
}

################# RUN #################

# 调用 CE_Sampling 函数时传递这些参数
with printto(log=log_path_01, lst=lst_path_01) as logger:
    CE1_Resampled, CE1_Sample_Rate = CE_Sampling(config=config, logger=logger)

with printto(log=log_path_02, lst=lst_path_02) as logger:
    CE2_Recoded, profile_df = CE_EDA_Recode(
        indsn=CE1_Resampled, config=config, logger=logger
    )
