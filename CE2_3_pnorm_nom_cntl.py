
import pandas as pd
import numpy as np
import os
from scipy.stats import t as student_t

# Define pnom and nom_cntl functions with corrected tmp computation

def pnom(df: pd.DataFrame,
         var: str,
         typ: bool,
         dep_var: str,
         vars2: pd.DataFrame,
         path_output: str,
         prefix: str,
         missrate: float,
         concrate: float,
         valcnt: int,
         minbinnc: int,
         minbinnp: float,
         talpha: float,
         bonfer: bool,
         nom_method: str,
         binary_dv: bool,
         overall_avg: float,
         nobs: int):
    # Summarize with explicit size, mean, var
    group = df.groupby(var, dropna=False)
    group_sizes = group.size().rename('dcount')
    group_means = group[dep_var].mean().rename('dmean')
    group_vars = group[dep_var].var().rename('dvar')
    tmp = pd.concat([group_sizes, group_means, group_vars], axis=1).reset_index().sort_values('dmean')
    # Counts
    misscnt = tmp.loc[tmp[var].isna(), 'dcount'].sum()
    maxcnt = tmp['dcount'].max()
    unqcnt = tmp.shape[0]
    # Update vars2
    cond = vars2['name'] == var
    vars2.loc[cond, ['unique_cnt','miss_cnt','max_cnt']] = [unqcnt, misscnt, maxcnt]
    # Checks
    if misscnt > nobs * missrate or maxcnt > nobs * concrate or (valcnt != 0 and unqcnt > valcnt):
        vars2.loc[cond, 'good'] = 0
        return vars2, pd.DataFrame(), []
    vars2.loc[cond, 'good'] = 1
    # Collapse by count
    if minbinnp > 0:
        minbinn = int(nobs * minbinnp)
    else:
        minbinn = minbinnc
    groups = []
    cumcnt = 0
    tcount = None
    tmean = None
    tvar = None
    group_id = 1
    for _, row in tmp.iterrows():
        cnt, mean, varr = row['dcount'], row['dmean'], row['dvar']
        cumcnt += cnt
        if tcount is None:
            tcount, tmean, tvar = cnt, mean, varr
        else:
            if tcount <= minbinn or (cumcnt - cnt) >= (nobs - minbinn):
                new_count = tcount + cnt
                new_mean = (tcount * tmean + cnt * mean) / new_count
                new_var = ((tcount - 1)*tvar + (cnt - 1)*varr) / (new_count - 2) if new_count > 2 else tvar
                tcount, tmean, tvar = new_count, new_mean, new_var
            else:
                group_id += 1
                tcount, tmean, tvar = cnt, mean, varr
        groups.append({'value': row[var], 'tgroup': group_id})
    # Bonferroni
    ncomps = group_id - 1
    ftalpha = talpha / ncomps if bonfer and ncomps > 1 else talpha
    # Collapse by variance
    gdf = pd.DataFrame(groups).merge(tmp[[var, 'dcount', 'dmean','dvar']], left_on='value', right_on=var)
    final_assign = {}
    fg = 1
    # Initialize first final group
    first = gdf[gdf['tgroup']==1]
    fcount = first['dcount'].sum()
    fmean = first['dmean'].iloc[0]
    fvar = first['dvar'].iloc[0]
    for _, grp in first.iterrows():
        final_assign[grp['value']] = fg
    for gid in range(2, group_id+1):
        group_rows = gdf[gdf['tgroup']==gid]
        tcount = group_rows['dcount'].sum()
        tmean = group_rows['dmean'].iloc[0]
        tvar = group_rows['dvar'].iloc[0]
        dfree = fcount + tcount - 2
        pvar = ((fcount-1)*fvar + (tcount-1)*tvar) / dfree if dfree > 0 else fvar
        t_stat = (fmean - tmean) / np.sqrt(pvar * (1/fcount + 1/tcount))
        p_val = 1 - student_t.cdf(abs(t_stat), dfree)
        if p_val <= ftalpha:
            fg += 1
            fcount, fmean, fvar = tcount, tmean, tvar
        else:
            new_count = fcount + tcount
            new_mean = (fcount*fmean + tcount*tmean) / new_count
            new_var = ((fcount-1)*fvar + (tcount-1)*tvar) / (new_count-2) if new_count>2 else fvar
            fcount, fmean, fvar = new_count, new_mean, new_var
        for _, grp in group_rows.iterrows():
            final_assign[grp['value']] = fg
    last_group = fg
    # Generate recode script
    recode_py = os.path.join(path_output, 'CE2_Nominal_Var_Recode.py')
    os.makedirs(path_output, exist_ok=True)
    with open(recode_py, 'a') as f:
        f.write(f"# Nominal recode for {var}\n")
        if nom_method.upper() == 'BINARY':
            for g in range(1, last_group):
                vals = [v for v, grp in final_assign.items() if grp == g]
                f.write(f"df['{prefix}{var}_X{g}'] = df['{var}'].isin({vals}).astype(int)\n")
        else:
            mapping = {v: df.loc[df[var]==v, dep_var].mean() for v in final_assign}
            f.write(f"df['{prefix}{var}'] = df['{var}'].map({mapping}).fillna(0)\n")
    # Profiling DataFrame
    profile = [{'variable': var, 'category': 'Overall',
                'Average_DV': overall_avg, 'Count': nobs,
                'Percent': 1.0, 'index': 100.0, 'star': ''}]
    group_vals = {}
    for val, grp in final_assign.items():
        group_vals.setdefault(grp, []).append(val)
    for grp, vals in group_vals.items():
        cnt = df[df[var].isin(vals)][dep_var].count()
        avg = df[df[var].isin(vals)][dep_var].mean()
        pct = cnt / nobs
        idx = (avg / overall_avg) * 100
        star = '* (+)' if idx >= 110 else '  (+)' if idx > 100 else '* (-)' if idx <= 90 else '  (-)'
        category = ','.join(map(str, vals)) if grp != last_group else 'Else'
        profile.append({'variable': var, 'category': category,
                        'Average_DV': avg, 'Count': cnt,
                        'Percent': pct, 'index': idx, 'star': star})
    profile_df = pd.DataFrame(profile)
    new_vars = ([f"{prefix}{var}_X{g}" for g in range(1, last_group)]
                if nom_method.upper()=='BINARY' else [prefix+var])
    return vars2, profile_df, new_vars

def nom_cntl(df: pd.DataFrame,
             nom_vars: list,
             vars2: pd.DataFrame,
             path_output: str,
             prefix: str,
             dep_var: str,
             missrate: float,
             concrate: float,
             valcnt: int,
             minbinnc: int,
             minbinnp: float,
             talpha: float,
             bonfer: bool,
             nom_method: str,
             binary_dv: bool):
    recode_py = os.path.join(path_output, 'CE2_Nominal_Var_Recode.py')
    os.makedirs(path_output, exist_ok=True)
    with open(recode_py, 'w') as f:
        f.write("# -*- coding: utf-8 -*-\n# Auto-generated nominal recode\n\n")
    profiles = []
    keep_list = []
    nobs = len(df)
    overall_avg = df[dep_var].mean()
    for var in nom_vars:
        vars2, prof, new_vars = pnom(df, var, df[var].dtype==int,
                                     dep_var, vars2, path_output,
                                     prefix, missrate, concrate,
                                     valcnt, minbinnc, minbinnp,
                                     talpha, bonfer, nom_method,
                                     binary_dv, overall_avg, nobs)
        if not prof.empty:
            profiles.append(prof)
        keep_list.extend(new_vars)
    with open(recode_py, 'a') as f:
        f.write("\n# KEEP_LIST_N\n")
        f.write(f"KEEP_LIST_N = {keep_list}\n")
    profile_df = pd.concat(profiles, ignore_index=True) if profiles else pd.DataFrame()
    return vars2, profile_df, keep_list