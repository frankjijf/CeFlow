
import os
import pandas as pd
import numpy as np

# 导入你已有的模块
from CE2_1_pnum_prof123 import pnum
from CE2_1_pnum_prof123 import prof1, prof2, prof3

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