
import os
import numpy as np
import pandas as pd

from CE2_1_pnum_prof123 import pnum
from CE2_1_pnum_prof123 import prof1, prof2, prof3

# 变换前缀到 Python 表达式的映射
_py_map = {
    'SQ_': lambda x: f"{x}**2",
    'SR_': lambda x: f"np.sqrt(np.maximum({x}, 0))",
    'LN_': lambda x: f"np.log(np.maximum({x}, 1e-5))",
    'IV_': lambda x: f"-1/np.maximum({x}, 1e-5)",
    'EP_': lambda x: f"-np.exp(np.minimum(-{x}, 0))",
    '':    lambda x: x
}

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