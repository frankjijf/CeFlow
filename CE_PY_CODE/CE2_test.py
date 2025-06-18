
import os
import numpy as np
import pandas as pd

from CE2_1_pnum_prof123 import pnum, prof1, prof2, prof3
from CE2_2_pbin_bin_cntl import bin_cntl
from CE2_3_pnorm_nom_cntl import nom_cntl
from CE2_4_pord_ord_cntl import ord_cntl
from CE2_5_pcon_cont_cntl import cont_cntl


# Define output directory and file paths
output_dir = "D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/"
titanic_csv_path = os.path.join(output_dir, "Titanic-Dataset.csv")
resampled_csv_path = os.path.join(output_dir, "CE1_Resampled_Titanic.csv")

# Load Titanic dataset
df_titanic = pd.read_csv(titanic_csv_path)

# Run CE_Sampling on Titanic (no bootstrap)
resampled_df, sample_rate_df = CE_Sampling(
    df=df_titanic,
    dep_var="Survived",
    split_frac=0.7,
    binary_dv=True,
    bootstrap=False,
    seed=654321,
    path_output=output_dir,
)

resampled_df.to_csv(resampled_csv_path, index=False)


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


# 假设 titanic 已经有 mod_val_test 列了
# 示例数据路径为硬编码，如需复现请将数据集放在相对路径下或修改为自己的数据路径
# 例如: df = pd.read_csv('./CE1_Resampled_Titanic.csv')
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/CE1_Resampled_Titanic.csv')

df['IsFemale'] = (df['Sex']=='female').astype(int)
df['IsChild'] = (df['Age'] < 18).astype(int)

# 定义各类变量列表
bin_vars  = ['IsFemale','IsChild']
nom_vars  = ['Embarked']
ord_vars  = ['Pclass']
cont_vars = ['Age','Fare']

recoded_df, profile_df, var_lookup_df = CE_EDA_Recode(
    insdn=df,
    bin_vars=bin_vars,
    nom_vars=nom_vars,
    ord_vars=ord_vars,
    cont_vars=cont_vars,
    prefix='R1_',
    dep_var='Survived',
    minbinnc=500,
    minbinnp=0.05,
    profiling=True,
    path_output='D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/'
)