
import os
import numpy as np
import pandas as pd

from CE2_1_pnum_prof123 import pnum, prof1, prof2, prof3
from CE2_2_pbin_bin_cntl import bin_cntl
from CE2_3_pnorm_nom_cntl import nom_cntl
from CE2_4_pord_ord_cntl import ord_cntl
from CE2_5_pcon_cont_cntl import cont_cntl

# Example usage of the pnum, prof1, prof2, and prof3 functions
# Load Titanic dataset
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/Titanic-Dataset.csv')

# Prepare vars2 DataFrame
vars2 = pd.DataFrame({
    'name': ['Age'],
    'var_lb': [np.nan],
    'var_ub': [np.nan],
    'var_median': [np.nan],
    'miss_impute': [np.nan],
    'Chisq': [np.nan],
    'PValue': [np.nan],
    'new_var': [''],
    'Sign': [''],
    'Relationship': ['']
})

# Test pnum on 'Age'
path_output = 'D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/'
vars2_updated = pnum(df=df,
                     var='Age',
                     dep_var='Survived',
                     typ='C',
                     vars2=vars2,
                     path_output=path_output,
                     prefix='R1_',
                     binary_dv=True,
                     transformationC=True,
                     cap_flrC=True)
# Display results
print(vars2_updated)
# Show recode script
with open(f"{path_output}/CE2_Continuous_Var_Recode.py", 'r') as f:
    script = f.read()

print("Generated SAS Recode Script:\n")
print(script)

# Calculate overall average and number of observations

overall_avg = df['Survived'].mean()
print(overall_avg)
nobs = len(df)
print(nobs)

# Test prof1 on 'Pclass'
prof1_pclass = prof1(df, 'Pclass', 'Survived')
print(prof1_pclass)

# Test prof2 on 'Age'
prof2_age = prof2(df, 'Age', 'Survived', num_category=5, equal_dist=False)
print(prof2_age)

# Prepare minimal vars2 for 'Age'
vars2_age = pd.DataFrame({'name': ['Age']})
print(vars2_age)

# Use prof2 output for prof3
prof3_age = prof3(vars2_age, prof2_age, var='Age', overall_avg=overall_avg, nobs=nobs, typ='C', dep_var='Survived')
print(prof3_age)

# Example usage
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/Titanic-Dataset.csv')

df['IsFemale'] = (df['Sex']=='female').astype(int)
df['IsChild'] = (df['Age'] < 18).astype(int)

# Prepare vars2
vars2 = pd.DataFrame({
    'name':['IsFemale','IsChild'],
    'label':['Is Female','Is Child']
})

# Test bin_cntl on both variables
vars2_updated, profile_df, keep_vars = bin_cntl(
    df, ['IsFemale','IsChild'], vars2.copy(),
    'D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/', prefix='R1_',
    dep_var='Survived', missrate=0.75,
    concrate=0.9, profiling=True,
    typ_map={'IsFemale':1,'IsChild':1}
)

# Display outputs
print(vars2_updated)
print(profile_df)
print(keep_vars)

# Print generated Python recode script
with open('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/CE2_Binary_Var_Recode.py','r') as f:
    script = f.read()
print("Generated CE2_Binary_Var_Recode.py:\n", script)


# Example usage:
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/Titanic-Dataset.csv')

# Prepare vars2 for 'Embarked'
vars2_nom = pd.DataFrame({
    'name': ['Embarked'],
    'label': ['Port of Embarkation']
})

# Test nom_cntl on 'Embarked'
vars2_updated, profile_nom, keep_nom = nom_cntl(
    df, ['Embarked'], vars2_nom.copy(),
    'D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/', 'R1_',
    'Survived', 0.75, 0.9, 50, 500, 0.05,
    0.05, True, 'BINARY', True
)

# Display results
print(vars2_updated)
print(profile_nom)
print(keep_nom)

# Show generated Python recode script
with open('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/CE2_Nominal_Var_Recode.py','r') as f:
    script = f.read()
print("Generated CE2_Nominal_Var_Recode.py:\n", script)


# Example usage of ord_cntl
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/Titanic-Dataset.csv')
vars2 = pd.DataFrame({
    'name': ['Age','Pclass'],
    'var_lb': [None,None],
    'var_ub': [None,None],
    'var_median': [None,None],
    'miss_impute': [None,None],
    'new_var': ['',''],
    'Sign': ['',''],
    'Relationship': ['',''],
    'unique_cnt': [None,None],
    'max_cnt': [None,None]
})
vars2_out, profile_df, keep_o = ord_cntl(
    df, ['Age','Pclass'], vars2,
    'D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/', 'R1_','Survived',
    missrate=0.75, concrate=0.9,
    profiling=True, transformationO=True,
    num_category=10
)

print("Updated vars2:\n", vars2_out)
print("Profile:\n", profile_df)
print("KEEP_LIST_O:", keep_o)


# 示例用法
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/Titanic-Dataset.csv')
# 初始 vars2 包含 col 'name'
vars2 = pd.DataFrame({'name': ['Age','Fare']})

vars2_out, profile_df, keep_c = cont_cntl(
    df,
    cont_vars=['Age','Fare'],
    vars2=vars2,
    path_output='D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/',
    prefix='R1_',
    dep_var='Survived',
    missrate=0.75,
    profiling=True,
    transformationC=True,
    p_lo=1,
    p_hi=99
)

print(vars2_out)
print(profile_df)
print("KEEP_LIST_C:", keep_c)
# 查看生成的 Python 脚本：
# /mnt/data/CE2_Continuous_Var_Recode.py
