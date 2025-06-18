
import pandas as pd
import numpy as np
import statsmodels.api as sm
import os
from scipy.stats import t as student_t
import glob
import logging
from typing import Optional

def pnum(df: pd.DataFrame,
         var: str,
         dep_var: str,
         typ: str,
         vars2: pd.DataFrame,
         path_output: str,
         prefix: str = '',
         missrate: float = 0.75,
         impmethodC: str = 'MEDIAN',
         impmethodO: str = 'MEAN',
         stdmethodC: str = 'STD',
         stdmethodO: str = 'NO',
         cap_flrC: bool = True,
         cap_flrO: bool = False,
         transformationC: bool = True,
         transformationO: bool = False,
         binary_dv: bool = True,
         min_size: int = 500,
         alpha: float = 0.05):

    # 1. Summary
    s = df[var]
    q = s.quantile([0.01,0.25,0.5,0.75,0.99])
    var_mean = s.mean()
    var_median = q[0.5]
    var_min, var_max = s.min(), s.max()
    var_p1, var_p25, var_p75, var_p99 = q[0.01], q[0.25], q[0.75], q[0.99]

    # 2. Bounds
    iqr = max(var_p75-var_p25, var_p99-var_p75, var_p25-var_p1)
    var_lb = min(max(var_p25-1.5*iqr, var_min), var_p1)
    var_ub = max(min(var_p75+1.5*iqr, var_max), var_p99)
    if var_lb==var_ub:
        var_lb, var_ub = var_min, var_max
    var_mid = (var_max-var_min)/2

    # 3. Default missing impute
    method = impmethodC.upper() if typ.upper()=='C' else impmethodO.upper()
    if method in ('MEAN','STD'):
        var_miss = var_mean
    elif method in ('MEDIAN','IQR','MAD'):
        var_miss = var_median
    elif method=='RANGE':
        var_miss = var_min
    elif method=='MIDRANGE':
        var_miss = var_mid
    elif method in ('SUM','EUCLEN','USTD','MAXABS'):
        var_miss = 0
    else:
        var_miss = None  # will compute for ER

    # 4. Cap/floor
    df_nm = df[[dep_var,var]].dropna()
    df_nm[var] = df_nm[var].clip(var_lb, var_ub)

    # 5. Transforms
    transforms = {}
    if (typ.upper()=='C' and transformationC) or (typ.upper()=='O' and transformationO):
        transforms[f'SQ_{var}'] = df_nm[var]**2
        transforms[f'SR_{var}'] = np.sqrt(np.maximum(df_nm[var],0))
        transforms[f'LN_{var}'] = np.log(np.maximum(df_nm[var],1e-5))
        transforms[f'IV_{var}'] = -1/np.maximum(df_nm[var],1e-5)
        transforms[f'EP_{var}'] = -np.exp(np.minimum(-df_nm[var],0))
    X = pd.DataFrame({var: df_nm[var], **transforms}); X = sm.add_constant(X)
    y = df_nm[dep_var]

    # 6. ER method stdize (not needed here as SAS uses PROC STDIZE only for non-ER)
    
    # 7. Univariate regression
    stats_map={}
    for col in [var]+list(transforms.keys()):
        df_x = X[['const', col]]
        if binary_dv:
            res = sm.Logit(y, df_x).fit(disp=False)
            stat = res.llr; p = res.pvalues[col]; est = res.params[col]; inter=res.params['const']
            stats_map[col]={'Intercept':inter,'Estimate':est,'Prob':p,'Stat':stat}
        else:
            res = sm.OLS(y, df_x).fit(disp=False)
            stat = res.fvalue; p = res.pvalues[col]; est = res.params[col]; inter=res.params['const']
            stats_map[col]={'Intercept':inter,'Estimate':est,'Prob':p,'Stat':stat,'RSquare':res.rsquared}

    # 8. Best transform
    best_col, best = min(stats_map.items(), key=lambda x:x[1]['Prob'])
    sign = '(+)' if best['Estimate']>0 else '(-)'
    pref = '' if best_col==var else best_col[:3]

    # 9. ER impute
    if method=='ER':
        miss_rate = df[df[var].isna()][dep_var].mean()
        if binary_dv:
            rr = np.clip(miss_rate,1e-4,0.9999); lo = np.log(rr/(1-rr))
            est, inter = best['Estimate'], best['Intercept']
            if pref=='SQ_': var_miss=np.sqrt(max((lo-inter)/est,0))
            elif pref=='SR_': var_miss=((lo-inter)/est)**2
            elif pref=='LN_': var_miss=np.exp((lo-inter)/est)
            elif pref=='IV_': var_miss=-1/((lo-inter)/est)
            elif pref=='EP_': var_miss=-np.exp(min(0,-(lo-inter)/est))
            else: var_miss=(lo-inter)/est
        else:
            mr = miss_rate; est, inter=best['Estimate'],best['Intercept']
            if pref=='SQ_': var_miss=np.sqrt(max((mr-inter)/est,0))
            elif pref=='SR_': var_miss=((mr-inter)/est)**2
            elif pref=='LN_': var_miss=np.exp((mr-inter)/est)
            elif pref=='IV_': var_miss=-1/((mr-inter)/est)
            elif pref=='EP_': var_miss=-np.exp(min(0,-(mr-inter)/est))
            else: var_miss=(mr-inter)/est
    var_miss=min(max(var_miss,var_lb),var_ub)

    # 10. Update vars2
    cond=vars2['name']==var
    vars2.loc[cond,['var_lb','var_ub','var_median','miss_impute']]=[var_lb,var_ub,var_median,var_miss]
    if binary_dv:
        vars2.loc[cond,['Chisq','PValue']]=[best['Stat'],best['Prob']]
    else:
        vars2.loc[cond,['FValue','PValue']]=[best['Stat'],best['Prob']]
    if best['Prob']<=alpha:
        vars2.loc[cond,['new_var','Sign','Relationship']]=[best_col,sign,pref+sign]

    # 11. Write Python recode
    py_fname = f"CE2_{'Continuous' if typ.upper()=='C' else 'Ordinal'}_Var_Recode.py"
    os.makedirs(path_output, exist_ok=True)
    op = os.path.join(path_output, py_fname)

    # 映射不同变换到 Python 表达式
    py_map = {
        'SQ_': lambda p: f"{p}**2",
        'SR_': lambda p: f"np.sqrt(np.maximum({p}, 0))",
        'LN_': lambda p: f"np.log(np.maximum({p}, 1e-5))",
        'IV_': lambda p: f"-1/np.maximum({p}, 1e-5)",
        'EP_': lambda p: f"-np.exp(np.minimum(-{p}, 0))",
        '':   lambda p: p
    }
    target = f"df['{prefix}{var}']"
    source = f"df['{var}']"

    with open(op, 'a') as f:
        f.write(f"# --- Recode variable: {var} ---\n")
        # ① 缺失填补
        f.write(f"{target} = {source}.fillna({var_miss})\n")
        # ② Capping/Flooring
        if (typ.upper()=='C' and cap_flrC) or (typ.upper()=='O' and cap_flrO):
            f.write(f"{target} = {target}.clip(lower={var_lb}, upper={var_ub})\n")
        # ③ Label 属性（可选）
        f.write(f"{target}.attrs['label'] = '{var}: Recode {sign}'\n")
        # ④ 最佳变换
        if pref:
            transform_expr = py_map[pref](target)
            f.write(f"df['{best_col}'] = {transform_expr}\n")
            f.write(f"df['{best_col}'].attrs['label'] = '{var} {pref} {sign}'\n")
        f.write("\n")

    return vars2

def prof1(df: pd.DataFrame,
          var: str,
          dep_var: str) -> pd.DataFrame:

    prof = (df
            .groupby(var)[dep_var]
            .agg(xcount='size', xmean='mean')
            .reset_index())
    prof['xcategory'] = prof[var].astype(str)
    prof.loc[df[var].isna(), 'xcategory'] = 'Missing'
    prof = prof.drop(columns=[var])
    return prof


def prof2(df: pd.DataFrame,
          var: str,
          dep_var: str,
          num_category: int = 10,
          equal_dist: bool = False) -> pd.DataFrame:

    tmp = df[[dep_var, var]].copy()
    if equal_dist:
        lo = df[var].quantile(0.01)
        hi = df[var].quantile(0.99)
        rng = (hi - lo) / num_category
        def assign_bin(x):
            if pd.isna(x): return np.nan
            if x < lo: return 1
            if x >= hi: return num_category
            return int((x - lo) // rng) + 1
        tmp['bin'] = tmp[var].map(assign_bin)
    else:
        tmp['bin'] = pd.qcut(tmp[var], q=num_category, labels=False, duplicates='drop') + 1
    prof = (tmp
            .groupby('bin', dropna=False)
            .agg(lo=(var, 'min'),
                 hi=(var, 'max'),
                 xcount=(var, 'size'),
                 xmean=(dep_var, 'mean'))
            .reset_index())
    prof['xcategory'] = prof.apply(
        lambda row: 'Missing' if pd.isna(row['bin']) else 
                   ('Low to ' + str(row['hi']) if row['bin']==1 else 
                    ('%s to High' % row['lo'] if row['bin']==prof['bin'].max() else
                     f"{row['lo']} to {row['hi']}")), axis=1)
    return prof


def prof3(vars2: pd.DataFrame,
          prof: pd.DataFrame,
          var: str,
          dep_var: str,
          overall_avg: float,
          nobs: int,
          typ: str) -> pd.DataFrame:

    merged = prof.copy()
    merged['variable'] = var
    merged = merged[['variable', 'xcategory', 'xcount', 'xmean']]
    out = []
    # Overall
    out.append({
        'variable': var,
        'category': 'Overall',
        'Average_DV': overall_avg,
        'Count': nobs,
        'Percent': 1.0,
        'index': 100,
        'star': ''
    })
    # Each bin
    for _, row in merged.iterrows():
        avg = row['xmean']
        pct = row['xcount'] / nobs
        idx = (avg / overall_avg) * 100 if overall_avg else np.nan
        if idx >= 110: star = '* (+)'
        elif idx > 100: star = '  (+)'
        elif idx <= 90: star = '* (-)'
        else: star = '  (-)'
        out.append({
            'variable': var,
            'category': row['xcategory'],
            'Average_DV': avg,
            'Count': row['xcount'],
            'Percent': pct,
            'index': idx,
            'star': star
        })
    result = pd.DataFrame(out)
    return result

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

# 变换前缀到Python表达式映射
_py_map = {
    'SQ_': lambda x: f"{x}**2",
    'SR_': lambda x: f"np.sqrt(np.maximum({x}, 0))",
    'LN_': lambda x: f"np.log(np.maximum({x}, 1e-5))",
    'IV_': lambda x: f"-1/np.maximum({x}, 1e-5)",
    'EP_': lambda x: f"-np.exp(np.minimum(-{x}, 0))",
    '':    lambda x: x
}

def ord_cntl(df: pd.DataFrame,
             ord_vars: list,
             vars2: pd.DataFrame,
             path_output: str,
             prefix: str,
             dep_var: str,
             missrate: float = 0.75,
             concrate: float = 0.9,
             profiling: bool = True,
             transformationO: bool = False,
             num_category: int = 10):
    """
    对应 SAS %ord_cntl：
      1. 初始化 Python recode 文件头
      2. 按 ord_vars 循环：
         a. 计算 missing count, unique_cnt, max_cnt, cons_dep → 更新 vars2
         b. 如果有效 → 调用 pord 处理该变量 （internally calls pnum，写 fillna/clip/transform）
         c. 按 need 调用 prof1/prof2/prof3，收集 profile_df
      3. 写 KEEP_LIST_O
      4. 返回：更新后 vars2, 合并 profile_df, keep_list
    """
    # 1. 初始化 recode 文件
    recode_py = os.path.join(path_output, "CE2_Ordinal_Var_Recode.py")
    os.makedirs(path_output, exist_ok=True)
    with open(recode_py, 'w') as f:
        f.write("# -*- coding: utf-8 -*-\n")
        f.write("# Auto-generated ordinal recode\n\n")

    all_profiles = []
    keep_list = []
    nobs = len(df)
    overall_avg = df[dep_var].mean()

    # 确保 vars2 包含所需列
    for col in ("unique_cnt","max_cnt","good","new_var","Sign","Relationship","miss_cnt"):
        if col not in vars2.columns:
            vars2[col] = np.nan

    # 2. 循环处理每个序数变量
    for var in ord_vars:
        # 2a. missing, unique, max, cons_dep
        miss_cnt = int(df[var].isna().sum())
        cnts = df[var].value_counts(dropna=False)
        unique_cnt = cnts.shape[0]
        max_cnt = int(cnts.max())
        nonmiss = df.loc[df[var].notna(), dep_var]
        cons_dep = nonmiss.nunique() == 1

        m = vars2['name'] == var
        vars2.loc[m, 'miss_cnt'] = miss_cnt
        vars2.loc[m, 'unique_cnt'] = unique_cnt
        vars2.loc[m, 'max_cnt'] = max_cnt

        # 2b. 测试
        if miss_cnt > nobs * missrate or max_cnt > nobs * concrate or cons_dep:
            # too many missing / too concentrated / constant DV → skip
            vars2.loc[m, 'good'] = 0
            continue
        vars2.loc[m, 'good'] = 1

        # 2c. 调用 pnum → 计算边界、impute、最佳变换 & 写出初始 recode
        # pnum 会写出 fillna、transform 部分并更新 vars2 中 var_lb,var_ub,miss_impute,new_var,Relationship,Sign
        vars2 = pnum(df=df,
                     var=var,
                     dep_var=dep_var,
                     typ='O',
                     vars2=vars2,
                     path_output=path_output,
                     prefix=prefix,
                     missrate=missrate,
                     transformationO=transformationO)

        # 2d. 更新 vars2 中的边界、impute、变换等信息
        meta      = vars2.loc[vars2['name']==var].iloc[0]
        var_lb    = meta['var_lb']
        var_ub    = meta['var_ub']
        var_miss  = meta['miss_impute']
        new_var   = meta.get('new_var','')    # 可能是 '', 或 'LN_Age' 等
        sign      = meta.get('Sign','')       # '(+)' 或 '(-)' 或 ''

        # 先写 clip & label
        with open(recode_py, 'a') as f:
            f.write(f"# --- Final clip & label for {var} ---\n")
            tgt = f"df['{prefix}{var}']"
            f.write(f"{tgt} = {tgt}.clip(lower={var_lb}, upper={var_ub})\n")
            f.write(f"{tgt}.attrs['label'] = '{var}: Recode {sign}'\n")

            # 再写 transform：根据 new_var 判断 prefix
            if new_var:
                valid_prefs = ('SQ_','SR_','LN_','IV_','EP_')
                # 如果 new_var 以某个前缀开头，就取它；否则没有前缀
                if any(new_var.startswith(pf) for pf in valid_prefs):
                    pref = new_var[:3]
                else:
                    pref = ''
                # 映射表
                mapping = {
                    'SQ_': lambda x: f"{x}**2",
                    'SR_': lambda x: f"np.sqrt(np.maximum({x}, 0))",
                    'LN_': lambda x: f"np.log(np.maximum({x}, 1e-5))",
                    'IV_': lambda x: f"-1/np.maximum({x}, 1e-5)",
                    'EP_': lambda x: f"-np.exp(np.minimum(-{x}, 0))",
                    '':    lambda x: x
                }
                expr = mapping[pref](tgt)
                f.write(f"df['{new_var}'] = {expr}\n")
                # 用 pref 去掉末尾的 '_' 做 label，若 pref='' 则只写 new_var
                lbl = pref[:-1] if pref else new_var
                f.write(f"df['{new_var}'].attrs['label'] = '{var} {lbl} {sign}'\n")

            f.write("\n")

        # 2e. Profiling
        if profiling:
            if unique_cnt <= num_category:
                tmp_prof = prof1(df, var, dep_var)
            else:
                tmp_prof = prof2(df, var, dep_var,
                                 num_category=num_category,
                                 equal_dist=False)
            prof3_df = prof3(vars2, tmp_prof,
                             var, dep_var,
                             overall_avg, nobs, 'O')
            all_profiles.append(prof3_df)

        # 记录 new_var 用于 KEEP_LIST_O
        if isinstance(new_var, str) and new_var.strip():
            keep_list.append(new_var)

    # 3. 写 KEEP_LIST_O
    with open(recode_py, 'a') as f:
        f.write("\n# KEEP_LIST_O\n")
        f.write(f"KEEP_LIST_O = {keep_list}\n")

    # 4. 合并所有 profiling 表
    profile_all = pd.concat(all_profiles, ignore_index=True) if all_profiles else pd.DataFrame()

    return vars2, profile_all, keep_list


def pcont(df: pd.DataFrame,
          var: str,
          vars2: pd.DataFrame,
          path_output: str,
          prefix: str,
          dep_var: str,
          missrate: float = 0.75,
          profiling: bool = True,
          transformationC: bool = True,
          p_lo: int = 1,
          p_hi: int = 99):
    """
    对应 SAS %pcont：
    - 检查常数 DV；否则
    - 调用 pnum(... typ='C') → 写 fillna & transform
    - 从 vars2 取出 metadata → 写 Python recode clip/label/transform
    - 如果 profiling=Y，做 prof2 + prof3('C') 并返回 profile_df
    返回 (vars2, profile_df, [new_var])
    """
    nobs = len(df)
    # 常数 DV 检查
    nonmiss_dep = df.loc[df[var].notna(), dep_var]
    if nonmiss_dep.nunique() <= 1:
        print(f"{var} has constant dependent value; excluded")
        return vars2, pd.DataFrame(), []

    # 调用 pnum 进行连续变量的核心计算与初步 Python recode
    vars2 = pnum(df=df,
                 var=var,
                 dep_var=dep_var,
                 typ='C',
                 vars2=vars2,
                 path_output=path_output,
                 prefix=prefix,
                 missrate=missrate,
                 transformationC=transformationC)

    # 从 vars2 中获取 recode 所需字段
    m = vars2['name'] == var
    meta = vars2.loc[m].iloc[0]
    var_lb    = meta['var_lb']
    var_ub    = meta['var_ub']
    var_miss  = meta['miss_impute']
    new_var   = meta.get('new_var', '')
    rel       = str(meta.get('Relationship',''))
    sign      = str(meta.get('Sign',''))
    pref      = rel[:3]  # e.g. 'LN_','SQ_'

    # 追加 Python recode：clip & label & transform
    recode_py = os.path.join(path_output, "CE2_Continuous_Var_Recode.py")
    os.makedirs(path_output, exist_ok=True)
    with open(recode_py, 'a') as f:
        f.write(f"# --- Recode continuous variable: {var} ---\n")
        tgt = f"df['{prefix}{var}']"
        src = f"df['{var}']"
        # 1. 填缺
        f.write(f"{tgt} = {src}.fillna({var_miss})\n")
        # 2. Capping/Flooring
        f.write(f"{tgt} = {tgt}.clip(lower={var_lb}, upper={var_ub})\n")
        # 3. Label
        f.write(f"{tgt}.attrs['label'] = '{var}: Recode {sign}'\n")
        # 4. 最佳变换
        if isinstance(new_var, str) and new_var:
            expr = _py_map[pref](tgt)
            f.write(f"df['{new_var}'] = {expr}\n")
            lbl = pref[:-1]
            f.write(f"df['{new_var}'].attrs['label'] = '{var} {lbl} {sign}'\n")
        f.write("\n")

    # Profiling: prof2 + prof3('C')
    profile_df = pd.DataFrame()
    if profiling:
        # prof2 分箱
        tmp = prof2(df, var, dep_var,
                    num_category=10,    # 用全局 num_category 或可传参
                    equal_dist=False)
        # prof3 合并
        overall_avg = df[dep_var].mean()
        profile_df = prof3(vars2, tmp,
                           var, dep_var,
                           overall_avg, nobs, 'C')

    return vars2, profile_df, [new_var] if new_var else []


def cont_cntl(df: pd.DataFrame,
              cont_vars: list,
              vars2: pd.DataFrame,
              path_output: str,
              prefix: str,
              dep_var: str,
              missrate: float = 0.75,
              profiling: bool = True,
              transformationC: bool = True,
              p_lo: int = 1,
              p_hi: int = 99):
    """
    对应 SAS %cont_cntl：
    - 初始化 CE2_Continuous_Var_Recode.py
    - 获取 vars2 基础 missing 等信息
    - 按 cont_vars 循环：miss / P_lo / P_hi / 检查 / 调 pcont
    - 写 KEEP_LIST_C
    返回 (vars2, combined_profile_df, keep_list)
    """
    # 初始化 recode 脚本
    recode_py = os.path.join(path_output, "CE2_Continuous_Var_Recode.py")
    os.makedirs(path_output, exist_ok=True)
    with open(recode_py,'w') as f:
        f.write("# -*- coding: utf-8 -*-\n")
        f.write("# Auto-generated continuous recode\n\n")

    all_profiles = []
    keep_list = []
    nobs = len(df)

    # 确保 vars2 包含 miss_cnt,P_lo,P_hi,good,new_var,Sign,Relationship
    for col in ("miss_cnt","P_lo","P_hi","good","new_var","Sign","Relationship"):
        if col not in vars2.columns:
            vars2[col] = np.nan

    # 循环
    for var in cont_vars:
        # 1. missing count
        miss_cnt = int(df[var].isna().sum())
        vars2.loc[vars2['name']==var, 'miss_cnt'] = miss_cnt

        # 2. P_lo / P_hi
        lo = df[var].quantile(p_lo/100)
        hi = df[var].quantile(p_hi/100)
        vars2.loc[vars2['name']==var, 'P_lo'] = lo
        vars2.loc[vars2['name']==var, 'P_hi'] = hi

        # 3. 测试：missing过多 or constant DV or lo==hi
        nonmiss_dep = df.loc[df[var].notna(), dep_var]
        if (miss_cnt > nobs * missrate) or (nonmiss_dep.nunique()<=1) or (lo==hi):
            vars2.loc[vars2['name']==var, 'good'] = 0
            continue
        vars2.loc[vars2['name']==var, 'good'] = 1

        # 4. 调用 pcont
        vars2, prof_df, new_vars = pcont(df, var, vars2,
                                         path_output, prefix,
                                         dep_var, missrate,
                                         profiling, transformationC,
                                         p_lo, p_hi)
        all_profiles.append(prof_df)
        keep_list += new_vars

    # 写 KEEP_LIST_C
    with open(recode_py,'a') as f:
        f.write("\n# KEEP_LIST_C\n")
        f.write(f"KEEP_LIST_C = {keep_list}\n")

    # 合并 profiling
    profile_all = pd.concat(all_profiles, ignore_index=True) if all_profiles else pd.DataFrame()

    return vars2, profile_all, keep_list


def CE_EDA_Recode(
    indsn: pd.DataFrame,
    config: dict,
    logger: Optional[logging.Logger] = None
):
    
    # Variable Lists
    binvar      = config["binvar"]
    nomvar      = config["nomvar"]    
    ordvar      = config["ordvar"] 
    contvar     = config["contvar"] 
    # General Macro Variables
    path_output   = config["path_output"] 
    inds          = config["inds"]
    id            = config["id"]
    dep_var       = config["dep_var"]
    binary_dv     = config["binary_dv"]
    weight     = config["weight"]
    # Optional Macro Variables: These all have defaults that can be used
    # General Macro Variables
    prefix = config.get("prefix", "R1_")
    keep_list = config.get("keep_list", [])
    # Macro 2: Recoding macro variables
    profiling = config.get("profiling", "Y").upper() == "Y"
    missrate = config.get("missrate", 0.75)
    concrate = config.get("concrate", 0.9)
    valcnt = config.get("valcnt", 50)
    minbinnc = config.get("minbinnc", 500)
    minbinnp = config.get("minbinnp", 0.05)
    talpha = config.get("talpha", 0.05)
    bonfer = config.get("bonfer", "N").upper() == "Y"
    nom_method = config.get("nom_method", "INDEX")
    pvalue = config.get("pvalue", 0.05)
    min_size = config.get("min_size", 500)
    num_category = config.get("num_category", 10)
    equal_dist = config.get("equal_dist", "N").upper() == "Y"
    p_lo = config.get("p_lo", 1)
    p_hi = config.get("p_hi", 99)
    impmethodC = config.get("impmethodC", "median")
    stdmethodC  = config.get("stdmethodC", "STD")
    cap_flrC = config.get("cap_flrC", "Y").upper() == "Y"
    transformationC = config.get("transformationC", "Y").upper() == "Y"
    impmethodO = config.get("impmethodO", "mean")
    stdmethodO  = config.get("stdmethodO", "No")
    cap_flrO = config.get("cap_flrO", "N").upper() == "Y"
    transformationO = config.get("transformationO", "N").upper() == "Y"

    # 1. Build workfile (exclude mod_val_test==3)
    workfile = pd.DataFrame(indsn.loc[indsn.get('mod_val_test', 0) != 3].copy())

    # 2. Global stats
    nobs = len(workfile)
    overall_avg = workfile[dep_var].mean()
    minbinn = max(minbinnc, int(nobs * minbinnp))

    # 3. Prepare profile DataFrame and var lookup
    profile_df = pd.DataFrame()
    var_lookup = []

    # 4. Binary variables
    if binvar:
        vars2_bin = pd.DataFrame({
            'name':  binvar,
            'label': binvar     # 这里给每个变量一个“原始标签”
        })
        typ_map = {v: (1 if pd.api.types.is_numeric_dtype(workfile[v]) else 2)
                for v in binvar}
        vars2_bin, prof_bin, keep_B = bin_cntl(
            df          = workfile,
            bin_vars    = binvar,
            vars2       = vars2_bin,
            path_output = path_output,
            prefix      = prefix,
            dep_var     = dep_var,
            missrate    = missrate,
            concrate    = concrate,
            profiling   = profiling,
            typ_map     = typ_map
        )
        if profiling and not prof_bin.empty:
            profile_df = pd.concat([profile_df, prof_bin], ignore_index=True)
    else:
        typ_map = {}
        vars2_bin = pd.DataFrame()
        keep_B = []

    # 5. Nominal variables
    if nomvar:
        vars2_nom = pd.DataFrame({
            'name':  nomvar,
            'label': nomvar
        })
        vars2_nom, prof_nom, keep_N = nom_cntl(
            df          = workfile,
            nom_vars    = nomvar,
            vars2       = vars2_nom,
            path_output = path_output,
            prefix      = prefix,
            dep_var     = dep_var,
            missrate    = missrate,
            concrate    = concrate,
            valcnt      = valcnt,
            minbinnc    = minbinnc,
            minbinnp    = minbinnp,
            talpha      = talpha,
            bonfer      = bonfer,
            nom_method  = nom_method,
            binary_dv   = binary_dv
        )
        if profiling and not prof_nom.empty:
            profile_df = pd.concat([profile_df, prof_nom], ignore_index=True)
    else:
        vars2_nom = pd.DataFrame()
        keep_N = []

    # 6. Ordinal variables
    if ordvar:
        vars2_ord = pd.DataFrame({
            'name':  ordvar,
            'label': ordvar
        })
        vars2_ord, prof_ord, keep_O = ord_cntl(
            df              = workfile,
            ord_vars        = ordvar,
            vars2           = vars2_ord,
            path_output     = path_output,
            prefix          = prefix,
            dep_var         = dep_var,
            missrate        = missrate,
            concrate        = concrate,
            profiling       = profiling,
            transformationO = transformationO,
            num_category    = num_category
        )
        if profiling and not prof_ord.empty:
            profile_df = pd.concat([profile_df, prof_ord], ignore_index=True)
    else:
        vars2_ord = pd.DataFrame()
        keep_O = []

    # 7. Continuous variables
    if contvar:
        vars2_cont = pd.DataFrame({
            'name':  contvar,
            'label': contvar
        })
        vars2_cont, prof_cont, keep_C = cont_cntl(
            df              = workfile,
            cont_vars       = contvar,
            vars2           = vars2_cont,
            path_output     = path_output,
            prefix          = prefix,
            dep_var         = dep_var,
            missrate        = missrate,
            profiling       = profiling,
            transformationC = transformationC,
            p_lo            = p_lo,
            p_hi            = p_hi
        )
        if profiling and not prof_cont.empty:
            profile_df = pd.concat([profile_df, prof_cont], ignore_index=True)
    else:
        vars2_cont = pd.DataFrame()
        keep_C = []

    # 8. Apply recode to entire dataset (ensuring new_vars exist)
    recoded = workfile.copy()
    for script in glob.glob(os.path.join(path_output, "CE2_*_Var_Recode.py")):
        # 把脚本读进来
        code = open(script, 'r').read()
        # 在 recoded DataFrame 上执行
        # 保证 df 名称一致，并且 np 也可用
        local_vars = {'df': recoded, 'np': np, 'nan': np.nan}
        exec(code, local_vars)
        # 更新 recoded，确保新变量被加入
        recoded = local_vars['df']
        exec(code, {'df': recoded, 'np': np, 'nan': np.nan})

    # 然后就能选出所有新变量
    base_keep = [c for c in [id, 'mod_val_test', dep_var] if c in recoded.columns]
    all_new   = keep_B + keep_N + keep_O + keep_C
    CE2_Recoded = recoded[base_keep + all_new].copy()

    # 10. Correlation and Excel report
    if dep_var in CE2_Recoded.columns and pd.api.types.is_numeric_dtype(CE2_Recoded[dep_var]):
        corr = CE2_Recoded.corr()[dep_var].drop(dep_var).abs().sort_values(ascending=False)
        corr_df = corr.reset_index().rename(columns={'index': 'variable', dep_var: 'correlation'})
    else:
        corr_df = pd.DataFrame(columns=['variable', 'correlation'])

    writer = pd.ExcelWriter(os.path.join(path_output, 'CE2_EDA_report.xlsx'),
                            engine='xlsxwriter')
    if keep_B: pd.DataFrame(vars2_bin).assign(new_var=pd.Series(keep_B)).to_excel(writer, 'Binary', index=False)
    if keep_N: pd.DataFrame(vars2_nom).assign(new_var=pd.Series(keep_N)).to_excel(writer, 'Nominal', index=False)
    if keep_O: pd.DataFrame(vars2_ord).assign(new_var=pd.Series(keep_O)).to_excel(writer, 'Ordinal', index=False)
    if keep_C: pd.DataFrame(vars2_cont).assign(new_var=pd.Series(keep_C)).to_excel(writer, 'Continuous', index=False)
    writer.close()

    # 11. Variable lookup
    lookup = []
    for v in all_new:
        lbl = v
        if not vars2_bin.empty and v in getattr(vars2_bin, 'name', []):
            try:
                lbl = vars2_bin.loc[vars2_bin['name'] == v, 'label'].values[0]
            except Exception:
                lbl = v
        elif not vars2_nom.empty and v in getattr(vars2_nom, 'name', []):
            try:
                lbl = vars2_nom.loc[vars2_nom['name'] == v, 'label'].values[0]
            except Exception:
                lbl = v
        elif not vars2_ord.empty and v in getattr(vars2_ord, 'name', []):
            try:
                lbl = vars2_ord.loc[vars2_ord['name'] == v, 'label'].values[0]
            except Exception:
                lbl = v
        elif not vars2_cont.empty and hasattr(vars2_cont, 'name') and v in vars2_cont['name'].values:
            try:
                lbl = vars2_cont.loc[vars2_cont['name'] == v, 'label'].values[0]
            except Exception:
                lbl = v
        lookup.append((v, lbl, v, lbl))
    var_lookup_df = pd.DataFrame(lookup,
                                 columns=['variable', 'label', 'orig_var', 'orig_label'])

    return CE2_Recoded, profile_df, var_lookup_df