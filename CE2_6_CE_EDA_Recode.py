
import os
import numpy as np
import pandas as pd

from CE2_2_pbin_bin_cntl import bin_cntl
from CE2_3_pnorm_nom_cntl import nom_cntl
from CE2_4_pord_ord_cntl import ord_cntl
from CE2_5_pcon_cont_cntl import cont_cntl

def CE_EDA_Recode(insdn: pd.DataFrame,
                  bin_vars: list,
                  nom_vars: list,
                  ord_vars: list,
                  cont_vars: list,
                  prefix: str,
                  dep_var: str,
                  minbinnc: int,
                  minbinnp: float,
                  profiling: bool,
                  path_output: str,
                  # 以下为各子宏默认参数
                  missrate: float = 0.75,
                  concrate: float = 0.9,
                  valcnt: int = 50,
                  talpha: float = 0.05,
                  bonfer: bool = False,
                  nom_method: str = 'INDEX',
                  transformationO: bool = False,
                  transformationC: bool = True,
                  num_category: int = 10,
                  p_lo: int = 1,
                  p_hi: int = 99):
    """
    Python 版 %CE_EDA_Recode，完整调用 bin_cntl/nom_cntl/ord_cntl/cont_cntl，
    并输出 recode 脚本与 EDA 报表，返回 recoded_df, profile_df, var_lookup_df。
    """
    # 1. Build workfile (exclude mod_val_test==3)
    workfile = insdn.loc[insdn.get('mod_val_test', 0) != 3].copy()

    # 2. Global stats
    nobs = len(workfile)
    overall_avg = workfile[dep_var].mean()
    minbinn = max(minbinnc, int(nobs * minbinnp))

    # 3. Prepare profile DataFrame and var lookup
    profile_df = pd.DataFrame()
    var_lookup = []

    # 4. Binary variables
    if bin_vars:
        vars2_bin = pd.DataFrame({
            'name':  bin_vars,
            'label': bin_vars     # 这里给每个变量一个“原始标签”
        })
        typ_map = {v: (1 if pd.api.types.is_numeric_dtype(workfile[v]) else 2)
                for v in bin_vars}
        vars2_bin, prof_bin, keep_B = bin_cntl(
            workfile, bin_vars, vars2_bin,
            path_output, prefix, dep_var,
            missrate, concrate, profiling, typ_map
        )
        if profiling and not prof_bin.empty:
            profile_df = pd.concat([profile_df, prof_bin], ignore_index=True)
    else:
        keep_B = []

    # 5. Nominal variables
    if nom_vars:
        vars2_nom = pd.DataFrame({
            'name':  nom_vars,
            'label': nom_vars
        })
        vars2_nom, prof_nom, keep_N = nom_cntl(
            workfile, nom_vars, vars2_nom,
            path_output, prefix, dep_var,
            missrate, concrate, valcnt, minbinnc, minbinnp,
            talpha, bonfer, nom_method,
            binary_dv=pd.api.types.is_numeric_dtype(insdn[dep_var]) and
                      insdn[dep_var].nunique() == 2
        )
        if profiling and not prof_nom.empty:
            profile_df = pd.concat([profile_df, prof_nom], ignore_index=True)
    else:
        keep_N = []

    # 6. Ordinal variables
    if ord_vars:
        vars2_ord = pd.DataFrame({
            'name':  ord_vars,
            'label': ord_vars
        })
        vars2_ord, prof_ord, keep_O = ord_cntl(
            workfile, ord_vars, vars2_ord,
            path_output, prefix, dep_var,
            missrate, concrate, profiling,
            transformationO, num_category
        )
        if profiling and not prof_ord.empty:
            profile_df = pd.concat([profile_df, prof_ord], ignore_index=True)
    else:
        keep_O = []

    # 7. Continuous variables
    if cont_vars:
        vars2_cont = pd.DataFrame({
            'name':  cont_vars,
            'label': cont_vars
        })
        vars2_cont, prof_cont, keep_C = cont_cntl(
            workfile, cont_vars, vars2_cont,
            path_output, prefix, dep_var,
            missrate, profiling, transformationC,
            p_lo, p_hi
        )
        if profiling and not prof_cont.empty:
            profile_df = pd.concat([profile_df, prof_cont], ignore_index=True)
    else:
        keep_C = []

    # 8. Apply recode to entire dataset (ensuring new_vars exist)
    recoded = insdn.copy()
    # run each cntl with profiling=False to generate new_vars on recoded
    if bin_vars:
        bin_cntl(
            df=recoded,
            bin_vars=bin_vars,
            vars2=vars2_bin,              # reuse the one with 'label'
            path_output=path_output,
            prefix=prefix,
            dep_var=dep_var,
            missrate=missrate,
            concrate=concrate,
            profiling=False,              # 这里关闭 profiling
            typ_map=typ_map               # 然后才是 typ_map
        )
    if nom_vars:
        nom_cntl(
            df=recoded,
            nom_vars=nom_vars,
            vars2=vars2_nom,
            path_output=path_output,
            prefix=prefix,
            dep_var=dep_var,
            missrate=missrate,
            concrate=concrate,
            valcnt=valcnt,
            minbinnc=minbinnc,
            minbinnp=minbinnp,
            talpha=talpha,
            bonfer=bonfer,
            nom_method=nom_method,
            binary_dv=(recoded[dep_var].nunique()==2)
        )
    if ord_vars:
        ord_cntl(
            df=recoded,
            ord_vars=ord_vars,
            vars2=vars2_ord,
            path_output=path_output,
            prefix=prefix,
            dep_var=dep_var,
            missrate=missrate,
            concrate=concrate,
            profiling=False,
            transformationO=transformationO,
            num_category=num_category
        )
    if cont_vars:
        cont_cntl(
            df=recoded,
            cont_vars=cont_vars,
            vars2=vars2_cont,
            path_output=path_output,
            prefix=prefix,
            dep_var=dep_var,
            missrate=missrate,
            profiling=False,
            transformationC=transformationC,
            p_lo=p_lo,
            p_hi=p_hi
        )

    # 9. Build CE2_Recoded (keep id, dep_var, mod_val_test and all new_vars)
    base_keep = [c for c in ['rid', 'mod_val_test', dep_var] if c in recoded.columns]
    all_new = keep_B + keep_N + keep_O + keep_C
    CE2_Recoded = recoded[base_keep + all_new].copy()

    # 10. Correlation and Excel report
    corr = CE2_Recoded.corr()[dep_var].drop(dep_var).abs().sort_values(ascending=False)
    corr_df = corr.reset_index().rename(columns={'index': 'variable', dep_var: 'correlation'})

    writer = pd.ExcelWriter(os.path.join(path_output, 'CE2_EDA_report.xlsx'),
                            engine='xlsxwriter')
    if keep_B: pd.DataFrame(vars2_bin).assign(new_var=keep_B).to_excel(writer, 'Binary', index=False)
    if keep_N: pd.DataFrame(vars2_nom).assign(new_vars=keep_N).to_excel(writer, 'Nominal', index=False)
    if keep_O: pd.DataFrame(vars2_ord).assign(new_var=keep_O).to_excel(writer, 'Ordinal', index=False)
    if keep_C: pd.DataFrame(vars2_cont).assign(new_var=keep_C).to_excel(writer, 'Continuous', index=False)
    if profiling and not profile_df.empty: profile_df.to_excel(writer, 'Profile', index=False)
    corr_df.to_excel(writer, 'Correlations', index=False)
    writer.save()

    # 11. Variable lookup
    lookup = []
    for v in all_new:
        lbl = recoded[v].attrs.get('label', '')
        lookup.append((v, lbl, v, lbl))
    var_lookup_df = pd.DataFrame(lookup,
                                 columns=['variable', 'label', 'orig_var', 'orig_label'])

    return CE2_Recoded, profile_df, var_lookup_df

# 假设 titanic 已经有 mod_val_test 列了
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