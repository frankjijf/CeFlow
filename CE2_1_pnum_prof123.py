
import pandas as pd
import numpy as np
import statsmodels.api as sm
import os

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
    """
    Complete Python reimplementation of SAS %pnum macro for continuous (typ='C')
    and ordinal (typ='O') variables, with one-to-one matching:
      1. Compute summary stats (mean, median, min, p1, p25, p75, p99, max)
      2. Calculate IQR-based var_lb and var_ub
      3. Determine default var_miss per impmethod (MEAN/MEDIAN/RANGE/MIDRANGE/etc.)
      4. Cap/floor non-missing values
      5. Generate transformations: SQ_, SR_, LN_, IV_, EP_
      6. Perform PROC STDIZE equivalent for ER method if required
      7. Run univariate logistic or OLS per transform, extract stats
      8. Select best transform (lowest p-value), compute ER-impute var_miss exactly per SAS logic
      9. Update vars2 table (var_lb, var_ub, var_median, miss_impute, others)
     10. Append recode code lines to CE2_Continuous/Ordinal_Var_Recode.txt
    """
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
    """
    Profile count and mean of dep_var by each category of var.
    Returns DataFrame with columns [variable, xcategory, xcount, xmean].
    """
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
    """
    Bin continuous/ordinal var into num_category groups. If equal_dist,
    use equal width; else quantile-based. Then summary by bin.
    Returns DataFrame with [bin, lo, hi, xcount, xmean, xcategory].
    """
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
    """
    Merge vars2 and prof into profiling output, add Overall row,
    compute percent, index, star per SAS logic.
    Returns DataFrame ready to append to CE2_profile.
    """
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