
import pandas as pd
import numpy as np
import os

# Function to process binary variables
# This function checks for missing values and counts of '1's in the binary variable,
# and generates a Python script to recode the variable.
def pbin(df: pd.DataFrame,
         var: str,
         typ: int,
         dep_var: str,
         vars2: pd.DataFrame,
         path_output: str,
         prefix: str,
         missrate: float,
         concrate: float,
         profiling: bool,
         overall_avg: float,
         nobs: int):
    miss_ck = df[var].isna().sum()
    miss_ok = (nobs * missrate) >= miss_ck
    cnt_ck = (df[var] == 1).sum() if typ == 1 else df[var].astype(str).str.strip().isin(['1','Y','y']).sum()
    lower, upper = nobs * (1 - concrate), nobs * concrate
    good = lower <= cnt_ck <= upper

    # Log the results
    recode_py = os.path.join(path_output, 'CE2_Binary_Var_Recode.py')
    with open(recode_py, 'a') as f:
        f.write(f"# Recode variable: {var}\n")
        f.write(f"df['{prefix}{var}'] = (df['{var}']==1).astype(int)\n")
        label = vars2.loc[vars2['name']==var,'label'].iloc[0]
        f.write(f"df['{prefix}{var}'].attrs['label'] = '{label}: Binary Recode'\n\n")

    # Profiling
    profile_df = None
    if profiling and miss_ok and good:
        tmp = df[[dep_var, var]].copy()
        tmp['newvar'] = (tmp[var] == 1).astype(int)
        prof = tmp.groupby('newvar')[dep_var].agg(xcount='size', xmean='mean').reset_index()
        rows = [{'variable': var, 'category': 'Overall',
                 'Average_DV': overall_avg, 'Count': nobs,
                 'Percent': 1.0, 'index': 100.0, 'star': ''}]
        for _, r in prof.iterrows():
            avg, cnt = r['xmean'], r['xcount']
            pct = cnt / nobs
            idx = (avg / overall_avg) * 100
            star = '* (+)' if idx >= 110 else '  (+)' if idx > 100 else '* (-)' if idx <= 90 else '  (-)'
            cat = '1' if r['newvar']==1 else 'Missing,0'
            rows.append({'variable': var, 'category': cat,
                         'Average_DV': avg, 'Count': cnt,
                         'Percent': pct, 'index': idx, 'star': star})
        profile_df = pd.DataFrame(rows)

    vars2.loc[vars2['name']==var, ['miss_cnt','ones_cnt','good']] = [miss_ck, cnt_ck, int(good and miss_ok)]
    return vars2, profile_df

def bin_cntl(df: pd.DataFrame,
             bin_vars: list,
             vars2: pd.DataFrame,
             path_output: str,
             prefix: str,
             dep_var: str,
             missrate: float,
             concrate: float,
             profiling: bool,
             typ_map: dict):
    # Prepare file
    os.makedirs(path_output, exist_ok=True)
    recode_py = os.path.join(path_output, 'CE2_Binary_Var_Recode.py')
    header = ["# -*- coding: utf-8 -*-",
              "# Auto-generated binary recode",
              ""]
    with open(recode_py, 'w') as f:
        f.write('\n'.join(header) + '\n')
    # Overall values
    nobs = len(df)
    overall_avg = df[dep_var].mean()
    profiles, new_vars = [], []
    # Loop
    for var in bin_vars:
        vars2, prof = pbin(df, var, typ_map[var], dep_var,
                           vars2, path_output, prefix,
                           missrate, concrate, profiling,
                           overall_avg, nobs)
        if prof is not None:
            profiles.append(prof)
        if vars2.loc[vars2['name']==var, 'good'].iloc[0] == 1:
            new_vars.append(prefix + var)
    # Write KEEP_LIST_B
    with open(recode_py, 'a') as f:
        f.write("\n# KEEP_LIST_B\n")
        f.write(f"KEEP_LIST_B = {new_vars}\n")
    profile_df = pd.concat(profiles, ignore_index=True) if profiles else pd.DataFrame()
    return vars2, profile_df, new_vars