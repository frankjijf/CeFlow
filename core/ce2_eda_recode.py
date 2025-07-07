import pandas as pd
import numpy as np
import statsmodels.api as sm
import os
from scipy.stats import t as student_t
import glob
import logging
from typing import Tuple, Optional
from sklearn.metrics import roc_auc_score

"""Utilities for exploratory data analysis and variable recoding.

This module contains helper functions used by the CEFlow pipeline to profile
variables, perform missing value imputation and create Python recode scripts.
Functions are grouped by variable type (binary, nominal, ordinal and
continuous) and return both updated metadata tables and optional profiling
DataFrames.
"""


def pnum(
    df: pd.DataFrame,
    var: str,
    typ: str,
    vars2: pd.DataFrame,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Profile a numeric variable and update the metadata table.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    var : str
        Variable name to process.
    typ : str
        Variable type ("C" or "O").
    vars2 : pd.DataFrame
        Table storing variable level information.
    config : dict
        Processing options.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    pd.DataFrame
        Updated ``vars2`` DataFrame.
    """

    dep_var = config["dep_var"]
    path_output = config["path_output"]
    binary_dv = config["binary_dv"].upper() == "Y"
    prefix = config.get("prefix", "R1_")
    impmethodC = config.get("impmethodC", "median")
    impmethodO = config.get("impmethodO", "mean")
    # stdmethodC = config.get('stdmethodC','STD')
    # stdmethodO = config.get('stdmethodO','No')
    # cap_flrC = config.get('cap_flrC','Y').upper() == 'Y'
    cap_flrO = config.get("cap_flrO", "N").upper() == "Y"
    transformationC = config.get("transformationC", "Y").upper() == "Y"
    transformationO = config.get("transformationO", "N").upper() == "Y"
    pvalue = config.get("pvalue", 0.05)
    min_size = config.get("min_size", 500)

    var_type = ""
    if typ.upper() == "O":
        var_type = "ordinal"
    if typ.upper() == "C":
        var_type = "continuous"
    if logger:
        logger.info(f"Processing numeric {var_type} variable: {var} in pnum")
        logger.debug(f"Processing numeric {var_type} variable: {var} in pnum")

    # 1. Summary
    s = df[var].dropna()
    q = s.quantile([0.01, 0.25, 0.5, 0.75, 0.99], interpolation="nearest")
    var_mean = s.mean()
    var_median = q[0.5]
    var_min = s.min()
    var_max = s.max()
    var_p1 = q[0.01]
    var_p25 = q[0.25]
    var_p75 = q[0.75]
    var_p99 = q[0.99]

    # 2. Bounds
    iqr = max(var_p75 - var_p25, var_p99 - var_p75, var_p25 - var_p1)
    var_lb = min(max(var_p25 - 1.5 * iqr, var_min), var_p1)
    var_ub = max(min(var_p75 + 1.5 * iqr, var_max), var_p99)
    if var_lb == var_ub:
        var_lb = var_min
        var_ub = var_max
    var_mid = (var_max - var_min) / 2

    # 3. Default missing impute
    skip = 1
    method = impmethodC.upper() if typ.upper() == "C" else impmethodO.upper()
    if method in ("MEAN", "STD"):
        var_miss = var_mean
    elif method in ("MEDIAN", "IQR", "MAD"):
        var_miss = var_median
    elif method == "RANGE":
        var_miss = var_min
    elif method == "MIDRANGE":
        var_miss = var_mid
    elif method in ("SUM", "EUCLEN", "USTD", "MAXABS"):
        var_miss = 0
    else:
        var_miss = None  # will compute for ER
        skip = 0

    # 4. Cap/floor
    df_nm = df[[dep_var, var]].dropna()
    df_nm.loc[:, var] = df_nm[var].clip(var_lb, var_ub)

    # 5. Transforms
    transforms = {}
    if (typ.upper() == "C" and transformationC) or (
        typ.upper() == "O" and transformationO
    ):
        transforms[f"SQ_{var}"] = df_nm[var] ** 2
        transforms[f"SR_{var}"] = np.sqrt(np.maximum(df_nm[var], 0))
        transforms[f"LN_{var}"] = np.log(np.maximum(df_nm[var], 0.00001))
        # transforms[f"IV_{var}"] = -1 / np.maximum(df_nm[var], 1e-5)
        # transforms[f"EP_{var}"] = -np.exp(np.minimum(-df_nm[var], 0))
    X = pd.DataFrame({var: df_nm[var], **transforms})
    if "const" not in X.columns:
        X = sm.add_constant(X)
    y = df_nm[dep_var]

    # 6. ER method stdize (in the future will update PROC STDIZE method for non-ER)
    # if skip = 0 THEN ABW, AHUBER, AWAVE, AGK, SPACING, L, and IN, TO BE UPDATED

    # 7. Univariate regression
    stats_map = {}
    for col in [var] + list(transforms.keys()):
        df_x = X[["const", col]]
        if binary_dv:
            res = sm.Logit(y, df_x).fit(disp=False)
            stat = res.llr
            p = res.pvalues[col]
            est = res.params[col]
            inter = res.params["const"]
            stats_map[col] = {
                "Intercept": inter,
                "Estimate": est,
                "Prob": p,
                "Stat": stat,
            }
        #            wald = res.wald_test_terms().statistic  # approximate WaldChiSq
        #            p_wald = res.pvalues[col]              # ProbChiSq
        #            auc = roc_auc_score(y, res.predict(df_x))
        #            stats_map[col].update({"WaldChiSq": wald, "ProbChiSq": p_wald, "CValue": auc})
        else:
            res = sm.OLS(y, df_x).fit(disp=False)
            stat = res.fvalue
            p = res.pvalues[col]
            est = res.params[col]
            inter = res.params["const"]
            stats_map[col] = {
                "Intercept": inter,
                "Estimate": est,
                "Prob": p,
                "Stat": stat,
                "RSquare": res.rsquared,
            }

    # 8. Best transform
    if stats_map:
        best_col, best = min(stats_map.items(), key=lambda x: x[1]["Prob"])
        sign = "(+)" if best["Estimate"] > 0 else "(-)"
        pref = "" if best_col == var else best_col[:3]
    else:
        best_col, best, sign, pref = (
            var,
            {"Estimate": 0, "Intercept": 0, "Prob": 1, "Stat": 0},
            "",
            "",
        )
    pref = "" if best_col == var else best_col[:3]

    # 9. For ER imputation ensure you use the dropna() subset; otherwise the denominator
    #    may differ when dep_var also has missing values.
    nmiss = df[var].isna().sum()
    ck = nmiss - min_size < 0
    ckP = best["Prob"] > 0.05
    if ck or ckP:
        var_miss = var_median
    else:
        if method == "ER":
            mrr_df = df.loc[df[var].isna() & df[dep_var].notna(), dep_var]
            miss_rate = (
                np.clip(mrr_df.mean(), 1e-4, 0.9999) if binary_dv else mrr_df.mean()
            )
            if binary_dv:
                rr = np.clip(miss_rate, 1e-4, 0.9999)
                lo = np.log(rr / (1 - rr))
                est, inter = best["Estimate"], best["Intercept"]
                if pref == "SQ_":
                    var_miss = np.sqrt(max((lo - inter) / est, 0))
                elif pref == "SR_":
                    var_miss = ((lo - inter) / est) ** 2
                elif pref == "LN_":
                    var_miss = np.exp((lo - inter) / est)
                elif pref == "IV_":
                    var_miss = -1 / ((lo - inter) / est)
                elif pref == "EP_":
                    var_miss = -np.exp(min(0, -(lo - inter) / est))
                else:
                    var_miss = (lo - inter) / est
            else:
                mr = miss_rate
                est, inter = best["Estimate"], best["Intercept"]
                if pref == "SQ_":
                    var_miss = np.sqrt(max((mr - inter) / est, 0))
                elif pref == "SR_":
                    var_miss = ((mr - inter) / est) ** 2
                elif pref == "LN_":
                    var_miss = np.exp((mr - inter) / est)
                elif pref == "IV_":
                    var_miss = -1 / ((mr - inter) / est)
                elif pref == "EP_":
                    var_miss = -np.exp(min(0, -(mr - inter) / est))
                else:
                    var_miss = (mr - inter) / est
        if var_miss is not None:
            var_miss = min(max(var_miss, var_lb), var_ub)
        else:
            # Fallback: assign a safe default if variables are uninitialized
            var_miss = var_median
        var_miss = min(max(var_miss, var_lb), var_ub)

    # 10. Update vars2
    target_cols = ["Chisq", "PValue", "FValue", "new_var", "Sign", "Relationship"]

    for col in target_cols:
        if col in vars2.columns:
            # If the column exists only change dtype without touching values
            vars2[col] = vars2[col].astype("string")
        else:
            # If missing create a full column of pd.NA with StringDtype
            vars2[col] = pd.Series(pd.NA, index=vars2.index, dtype="string")

    cond = vars2["var_name"] == var
    vars2.loc[cond, "var_lb"] = var_lb
    vars2.loc[cond, "var_ub"] = var_ub
    vars2.loc[cond, "var_median"] = var_median
    vars2.loc[cond, "miss_impute"] = var_miss

    if binary_dv:
        vars2.loc[cond, "Chisq"] = str(best["Stat"])
        vars2.loc[cond, "PValue"] = str(best["Prob"])
    else:
        vars2.loc[cond, "FValue"] = str(best["Stat"])
        vars2.loc[cond, "PValue"] = str(best["Prob"])
    if best["Prob"] <= pvalue:
        vars2.loc[cond, "new_var"] = str(best_col)
        vars2.loc[cond, "Sign"] = str(sign)
        vars2.loc[cond, "Relationship"] = str(pref + sign)

    # 11. Write Python recode
    os.makedirs(path_output, exist_ok=True)
    py_fname = f"CE2_{'Continuous' if typ.upper()=='C' else 'Ordinal'}_Var_Recode.py"
    recode_py = os.path.join(path_output, py_fname)

    # 2) Open file and begin writing
    with open(recode_py, "a", encoding="utf-8") as f:
        f.write(f"# --- Recode variable: {var} ---\n")
        # 3) Missing imputation: explicit IF/ELSE
        #    For missing values
        f.write(f"df.loc[df['{var}'].isna(), '{prefix}{var}'] = {var_miss}\n")
        #    For non-missing values
        f.write(f"df.loc[df['{var}'].notna(), '{prefix}{var}'] = df['{var}']\n")
        # 4) Capping/Flooring
        if cap_flrO:
            f.write(
                f"df['{prefix}{var}'] = "
                f"df['{prefix}{var}'].clip(lower={var_lb}, upper={var_ub})\n"
            )
        # 5) Create new column using the best transform
        if pref:  # e.g. 'SQ_', 'LN_', etc.
            # New column name
            new_col = f"{prefix}{pref}{var}"
            # Expression
            expr = {
                "SQ_": f"df['{prefix}{var}']**2",
                "SR_": f"np.sqrt(np.maximum(df['{prefix}{var}'], 0))",
                "LN_": f"np.log(np.maximum(df['{prefix}{var}'], 1e-5))",
                "IV_": f"-1/np.maximum(df['{prefix}{var}'], 1e-5)",
                "EP_": f"-np.exp(np.minimum(-df['{prefix}{var}'], 0))",
            }[pref]
            f.write(f"df['{new_col}'] = {expr}\n")
        f.write("\n")
        if logger:
            logger.info(f"Recode for variable {var} written to {recode_py}")
            logger.debug(f"Recode for variable {var} written to {recode_py}")

    return vars2


def prof1(
    df: pd.DataFrame,
    var: str,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Profile ``var`` by calculating count and mean of the dependent variable.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset containing ``var`` and the dependent variable.
    var : str
        Column name to group by.
    config : dict
        Configuration dictionary with ``dep_var`` specifying the dependent
        variable.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ``xcategory``, ``xcount`` and ``xmean``.
    """

    dep_var = config["dep_var"]
    prof = df.groupby(var)[dep_var].agg(xcount="size", xmean="mean").reset_index()
    prof["xcategory"] = prof[var].astype(str)
    prof.loc[df[var].isna(), "xcategory"] = "Missing"
    prof = prof.drop(columns=[var])
    return prof


def prof2(
    df: pd.DataFrame,
    var: str,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Bin a numeric variable and profile the dependent variable in each bin.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    var : str
        Continuous variable to bin.
    config : dict
        Configuration dictionary with ``dep_var`` and optional binning options
        like ``num_category`` and ``equal_dist``.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    pd.DataFrame
        Summary statistics for each bin including ``lo``, ``hi``, ``xcount`` and
        ``xmean``.
    """

    dep_var = config["dep_var"]
    num_category = config.get("num_category", 10)
    equal_dist = config.get("equal_dist", "N") == "Y"
    tmp = df[[dep_var, var]].copy()
    if equal_dist:
        lo = df[var].quantile(0.01)
        hi = df[var].quantile(0.99)
        rng = (hi - lo) / num_category

        def assign_bin(x):
            if pd.isna(x):
                return np.nan
            if x < lo:
                return 1
            if x >= hi:
                return num_category
            return int((x - lo) // rng) + 1

        tmp["bin"] = tmp[var].map(assign_bin)
    else:
        tmp["bin"] = (
            pd.qcut(tmp[var], q=num_category, labels=False, duplicates="drop") + 1
        )
    prof = (
        tmp.groupby("bin", dropna=False)
        .agg(
            lo=(var, "min"),
            hi=(var, "max"),
            xcount=(var, "size"),
            xmean=(dep_var, "mean"),
        )
        .reset_index()
    )
    prof["xcategory"] = prof.apply(
        lambda row: (
            "Missing"
            if pd.isna(row["bin"])
            else (
                "Low to " + str(row["hi"])
                if row["bin"] == 1
                else (
                    "%s to High" % row["lo"]
                    if row["bin"] == prof["bin"].max()
                    else f"{row['lo']} to {row['hi']}"
                )
            )
        ),
        axis=1,
    )
    return prof


def prof3(
    prof: pd.DataFrame,
    var: str,
    overall_avg: float,
    nobs: int,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """Convert profiling results to a standardized report format.

    Parameters
    ----------
    prof : pd.DataFrame
        Output from :func:`prof1` or :func:`prof2`.
    var : str
        Variable name that was profiled.
    overall_avg : float
        Overall mean of the dependent variable.
    nobs : int
        Total number of observations used in profiling.
    config : dict
        Configuration dictionary (currently unused).
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    pd.DataFrame
        DataFrame summarising each category/bin with percentage, index and
        highlight star.
    """

    merged = prof.copy()
    merged["variable"] = var
    merged = merged[["variable", "xcategory", "xcount", "xmean"]]
    out = []
    # Overall
    out.append(
        {
            "Variable": var,
            "Category": "Overall",
            "Average_DV": overall_avg,
            "Count": nobs,
            "Percent": 1.0,
            "Index": 100,
            "Star": "",
        }
    )
    # Each bin
    for _, row in merged.iterrows():
        avg = row["xmean"]
        pct = row["xcount"] / nobs
        idx = (avg / overall_avg) * 100 if overall_avg else np.nan
        star = (
            "* (+)"
            if idx >= 110
            else (
                "  (+)"
                if idx > 100
                else "  (-)" if idx > 90 else "* (-)" if idx <= 90 else "  (0)"
            )
        )
        out.append(
            {
                "Variable": var,
                "Category": row["xcategory"],
                "Average_DV": avg,
                "Count": row["xcount"],
                "Percent": pct,
                "Index": idx,
                "Star": star,
            }
        )
    result = pd.DataFrame(out)
    return result


def pbin(
    df: pd.DataFrame,
    var: str,
    typ: str,
    vars2_bin: pd.DataFrame,
    overall_avg: float,
    nobs: int,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process a binary variable and generate recode instructions.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    var : str
        Variable to process.
    typ : str
        ``"numeric"`` if values are 0/1, otherwise ``"str"``.
    vars2_bin : pd.DataFrame
        Metadata table to update.
    overall_avg : float
        Overall average of the dependent variable.
    nobs : int
        Number of observations in ``df``.
    config : dict
        Configuration dictionary with options such as ``missrate`` and
        ``concrate``.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Updated ``vars2_bin`` and an optional profiling DataFrame.
    """

    dep_var = config["dep_var"]
    path_output = config["path_output"]
    prefix = config.get("prefix", "R1_")
    missrate = config.get("missrate", 0.75)
    concrate = config.get("concrate", 0.9)
    profiling = config.get("profiling", "Y") == "Y"

    if logger:
        logger.info(f"Processing binary variable: {var}, type: {typ} in pbin")
        logger.debug(f"Processing binary variable: {var}, type: {typ} in pbin")
    prof_bin = pd.DataFrame()
    miss_cnt = df[var].isna().sum()
    miss_ok = (nobs * missrate) >= miss_cnt
    if miss_ok:
        if typ == "numeric":
            cnt_ck = (df[var] == 1).sum()
        else:
            cnt_ck = (
                df[var]
                .astype(str)
                .str.strip()
                .isin(["True", "TRUE", "1", "Y", "y"])
                .sum()
            )
        concrate_upper = nobs * concrate
        concrate_lower = nobs * (1 - concrate)
        concrate_good = concrate_lower <= cnt_ck <= concrate_upper
        if concrate_good:
            if logger:
                logger.info(
                    f"Processing binary variable: {var}, type: {typ}, miss_ok: {miss_ok}, cnt_ck: {cnt_ck}, concrate_good: {concrate_good}"
                )
                logger.debug(
                    f"Processing binary variable: {var}, type: {typ}, miss_ok: {miss_ok}, cnt_ck: {cnt_ck}, concrate_good: {concrate_good}"
                )
            os.makedirs(path_output, exist_ok=True)
            recode_py = os.path.join(path_output, "CE2_Binary_Var_Recode.py")
            with open(recode_py, "a", encoding="utf-8") as f:
                f.write(f"# Recode variable: {var}\n")
                if typ == "numeric":
                    recode_expr = f"(df['{var}'] == 1).fillna(False).astype(int)"
                else:
                    truth_values = "['True','TRUE','1','Y','y']"
                    recode_expr = (
                        f"df['{var}']"
                        f".astype(str)"
                        f".str.strip()"
                        f".isin({truth_values})"
                        f".astype(int)"
                    )
                f.write(f"df['{prefix}{var}'] = {recode_expr}\n")
            if logger:
                logger.info(f"Recode for variable {var} written to {recode_py}")
                logger.debug(f"Recode for variable {var} written to {recode_py}")
            # update vars2_bin
            cond = vars2_bin["var_name"] == var
            vars2_bin.loc[cond, ["miss_cnt", "ones_cnt", "overall_good"]] = [
                miss_cnt,
                cnt_ck,
                int(miss_ok and concrate_good),
            ]
            if logger:
                logger.info(
                    f"Updated vars2_bin for {var}: miss_cnt={miss_cnt}, ones_cnt={cnt_ck}, overall_good={int(miss_ok and concrate_good)}"
                )
                logger.debug(
                    f"Updated vars2_bin for {var}: miss_cnt={miss_cnt}, ones_cnt={cnt_ck}, overall_good={int(miss_ok and concrate_good)}"
                )
            # Profiling
            if profiling:
                tmp = df[[dep_var, var]].copy()
                if typ == "numeric":
                    tmp["newvar"] = (tmp[var] == 1).fillna(False).astype(int)
                else:
                    truth_values = [
                        "True",
                        "TRUE",
                        "1",
                        "Y",
                        "y",
                        "Yes",
                        "yes",
                        "On",
                        "on",
                    ]
                    tmp["newvar"] = (
                        tmp[var]
                        .astype(str)
                        .str.strip()
                        .isin(truth_values)
                        .fillna(False)
                        .astype(int)
                    )
                prof = (
                    tmp.groupby("newvar")[dep_var]
                    .agg(xcount="size", xmean="mean")
                    .reset_index()
                )
                profiles = [
                    {
                        "Variable": var,
                        "Category": "Overall",
                        "Average_DV": overall_avg,
                        "Count": nobs,
                        "Percent": 1.0,
                        "Index": 100.0,
                        "Star": "",
                    }
                ]
                for _, r in prof.iterrows():
                    avg, cnt = r["xmean"], r["xcount"]
                    pct = cnt / nobs
                    idx = (avg / overall_avg) * 100
                    star = (
                        "* (+)"
                        if idx >= 110
                        else (
                            "  (+)"
                            if idx > 100
                            else (
                                "  (-)"
                                if idx > 90
                                else "* (-)" if idx <= 90 else "  (0)"
                            )
                        )
                    )
                    cat = "1" if r["newvar"] == 1 else "Missing,0"
                    profiles.append(
                        {
                            "Variable": var,
                            "Category": cat,
                            "Average_DV": avg,
                            "Count": cnt,
                            "Percent": pct,
                            "Index": idx,
                            "Star": star,
                        }
                    )
                prof_bin = pd.DataFrame(profiles)
                if logger:
                    logger.info(f"Profiling for variable {var} completed.")
                    logger.debug(f"Profiling for variable {var} completed.")
        else:
            # update vars2_bin
            cond = vars2_bin["var_name"] == var
            vars2_bin.loc[cond, ["miss_cnt", "ones_cnt", "overall_good"]] = [
                miss_cnt,
                cnt_ck,
                int(miss_ok and concrate_good),
            ]
            if logger:
                logger.warning(
                    f"Variable {var} is too concentrated to use: {cnt_ck} not in [{concrate_lower}, {concrate_upper}]"
                )
                logger.debug(
                    f"Variable {var} is too concentrated to use: {cnt_ck} not in [{concrate_lower}, {concrate_upper}]"
                )
    else:
        # update vars2_bin
        cond = vars2_bin["var_name"] == var
        vars2_bin.loc[cond, ["miss_cnt", "ones_cnt", "overall_good"]] = [
            miss_cnt,
            np.nan,
            0,
        ]
        if logger:
            logger.warning(
                f"Variable {var} missing count is {miss_cnt}, which is too high"
            )
            logger.debug(
                f"Variable {var} missing count is {miss_cnt}, which is too high"
            )
    return vars2_bin, prof_bin


def bin_cntl(
    df: pd.DataFrame,
    bin_vars: list,
    vars2_bin: pd.DataFrame,
    typ_map: dict,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Controller for processing multiple binary variables.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    bin_vars : list
        List of binary variable names to process.
    vars2_bin : pd.DataFrame
        Metadata table for binary variables.
    typ_map : dict
        Mapping of variable name to ``"numeric"`` or ``"str"``.
    config : dict
        Configuration dictionary used by :func:`pbin`.
    logger : logging.Logger, optional
        Optional logger for progress output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Updated ``vars2_bin``, profiling DataFrame and keep list for recoded
        variables.
    """

    path_output = config["path_output"]
    prefix = config.get("prefix", "R1_")
    dep_var = config["dep_var"]

    if logger:
        logger.info(f"bin_cntl start, number of binary variables = {len(bin_vars)}")
        logger.debug(f"bin_cntl start, number of binary variables = {len(bin_vars)}")
    #
    os.makedirs(path_output, exist_ok=True)
    recode_py = os.path.join(path_output, "CE2_Binary_Var_Recode.py")
    header = ["# -*- coding: utf-8 -*-", "# Auto-generated binary recode", ""]
    with open(recode_py, "w") as f:
        f.write("\n".join(header) + "\n")
    #
    prof_bin_all = []
    new_vars = []
    nobs = len(df)
    overall_avg = df[dep_var].mean()
    #
    for var in bin_vars:
        vars2_bin, prof_bin = pbin(
            df=df,
            var=var,
            typ=typ_map[var],
            vars2_bin=vars2_bin,
            overall_avg=overall_avg,
            nobs=nobs,
            config=config,
            logger=logger,
        )
        if prof_bin is not None:
            prof_bin_all.append(prof_bin)
        vars2_bin.loc[vars2_bin["var_name"] == var, "new_var"] = ""
        if vars2_bin.loc[vars2_bin["var_name"] == var, "overall_good"].iloc[0] == 1:
            vars2_bin.loc[vars2_bin["var_name"] == var, "new_var"] = prefix + var
            new_vars.append(prefix + var)
        if logger:
            logger.info(
                f"Processed variable: {var}, overall_good={vars2_bin.loc[vars2_bin['var_name'] == var, 'overall_good'].iloc[0]}"
            )
            logger.debug(
                f"Processed variable: {var}, overall_good={vars2_bin.loc[vars2_bin['var_name'] == var, 'overall_good'].iloc[0]}"
            )
    #
    with open(recode_py, "a") as f:
        f.write("\n# KEEP_LIST_B\n")
        f.write(f"KEEP_LIST_B = {new_vars}\n")

    keep_vars_bin = vars2_bin[vars2_bin["overall_good"] == 1][
        ["var_name", "new_var"]
    ].copy()
    keep_vars_bin = keep_vars_bin.rename(
        columns={"var_name": "orig_var", "new_var": "variable"}
    )

    profile_df = (
        pd.concat(prof_bin_all, ignore_index=True) if prof_bin_all else pd.DataFrame()
    )

    if logger:
        logger.info(f"bin_cntl end, number of new binary variables = {len(new_vars)}")
        logger.debug(f"bin_cntl end, number of new binary variables = {len(new_vars)}")

    return vars2_bin, profile_df, keep_vars_bin


def pnom(
    df: pd.DataFrame,
    var: str,
    typ: str,
    dep_var: str,
    vars2_nom: pd.DataFrame,
    overall_avg: float,
    nobs: int,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process a nominal (categorical) variable.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    var : str
        Variable name to process.
    typ : str
        ``"numeric"`` or ``"str`` describing the data type of ``var``.
    dep_var : str
        Name of the dependent variable column.
    vars2_nom : pd.DataFrame
        Metadata table for nominal variables.
    overall_avg : float
        Overall average of the dependent variable.
    nobs : int
        Number of observations.
    config : dict
        Configuration dictionary controlling grouping and thresholds.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Updated ``vars2_nom`` and a profiling DataFrame for ``var``.
    """

    path_output = config["path_output"]
    prefix = config.get("prefix", "R1_")
    missrate = config.get("missrate", 0.75)
    concrate = config.get("concrate", 0.9)
    profiling = config.get("profiling", "Y") == "Y"
    valcnt = config.get("valcnt", 50)
    minbinnc = config.get("minbinnc", 500)
    minbinnp = config.get("minbinnp", 0.05)
    talpha = config.get("talpha", 0.05)
    bonfer = config.get("bonfer", "N") == "Y"
    nom_method = config.get("nom_method", "INDEX")

    if logger:
        logger.info(f"Processing nominal variable pnom: {var}, type: {typ}")
        logger.debug(f"Processing nominal variable pnom: {var}, type: {typ}")

    prof_nom = pd.DataFrame()
    # Summarize with explicit size, mean, var
    group = df.groupby(var, dropna=False)
    group_sizes = group.size().rename("dcount")
    group_means = group[dep_var].mean().rename("dmean")
    group_vars = group[dep_var].var().rename("dvar")
    tmp = (
        pd.concat([group_sizes, group_means, group_vars], axis=1)
        .reset_index()
        .sort_values("dmean")
    )
    #
    miss_cnt = tmp.loc[tmp[var].isna(), "dcount"].sum()
    miss_ok = (nobs * missrate) >= miss_cnt
    if not miss_ok:
        if logger:
            logger.warning(
                f"Variable {var} missing count is {miss_cnt}, which is too high"
            )
            logger.debug(
                f"Variable {var} missing count is {miss_cnt}, which is too high"
            )
        # Update vars2_nom
        cond = vars2_nom["var_name"] == var
        vars2_nom.loc[cond, ["uniq_cnt", "miss_cnt", "max_cnt"]] = [
            np.nan,
            miss_cnt,
            np.nan,
        ]
    else:
        max_cnt = tmp["dcount"].max()
        concrate_lower = nobs * (1 - concrate)
        concrate_upper = nobs * concrate
        concrate_good = concrate_lower <= max_cnt <= concrate_upper
        if not concrate_good:
            if logger:
                logger.warning(
                    f"Variable {var} is too concentrated to use: {max_cnt} not in [{concrate_lower}, {concrate_upper}]"
                )
                logger.debug(
                    f"Variable {var} is too concentrated to use: {max_cnt} not in [{concrate_lower}, {concrate_upper}]"
                )
            # Update vars2_nom
            cond = vars2_nom["var_name"] == var
            vars2_nom.loc[cond, ["uniq_cnt", "miss_cnt", "max_cnt"]] = [
                np.nan,
                miss_cnt,
                max_cnt,
            ]
        else:
            uniq_cnt = tmp.shape[0]
            val_ok = valcnt != 0 and uniq_cnt <= valcnt
            if not val_ok:
                if logger:
                    logger.warning(
                        f"Variable {var} has has too many values to use: {uniq_cnt} exceeds the limit of {valcnt}"
                    )
                    logger.debug(
                        f"Variable {var} has has too many values to use: {uniq_cnt} exceeds the limit of {valcnt}"
                    )
                # Update vars2_nom
                cond = vars2_nom["var_name"] == var
                vars2_nom.loc[cond, ["uniq_cnt", "miss_cnt", "max_cnt"]] = [
                    uniq_cnt,
                    miss_cnt,
                    max_cnt,
                ]
            else:
                # Update vars2_nom
                cond = vars2_nom["var_name"] == var
                vars2_nom.loc[cond, ["uniq_cnt", "miss_cnt", "max_cnt"]] = [
                    uniq_cnt,
                    miss_cnt,
                    max_cnt,
                ]
                cond = vars2_nom["var_name"] == var
                vars2_nom.loc[cond, ["overall_good"]] = 1

            # 1. Collapse values based on counts
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
                cnt, mean, varr = row["dcount"], row["dmean"], row["dvar"]
                cumcnt += cnt
                if tcount is None:
                    tcount, tmean, tvar = cnt, mean, varr
                else:
                    if tcount <= minbinn or (cumcnt - cnt) >= (nobs - minbinn):
                        new_count = tcount + cnt
                        new_mean = (tcount * tmean + cnt * mean) / new_count
                        new_var = (
                            ((tcount - 1) * tvar + (cnt - 1) * varr) / (new_count - 2)
                            if new_count > 2
                            else tvar
                        )
                        tcount, tmean, tvar = new_count, new_mean, new_var
                    else:
                        group_id += 1
                        tcount, tmean, tvar = cnt, mean, varr
                groups.append({"value": row[var], "tgroup": group_id})

            # 2. Calculate bonferroni adjustment on alpha value
            ncomps = group_id - 1
            ftalpha = talpha / ncomps if bonfer and ncomps > 1 else talpha

            # 3. Collapse values based on variance
            gdf = pd.DataFrame(groups).merge(
                tmp[[var, "dcount", "dmean", "dvar"]], left_on="value", right_on=var
            )
            final_assign = {}
            fg = 1
            # Initialize first final group
            first = gdf[gdf["tgroup"] == 1]
            fcount = first["dcount"].sum()
            fmean = first["dmean"].iloc[0]
            fvar = first["dvar"].iloc[0]
            for _, grp in first.iterrows():
                final_assign[grp["value"]] = fg
            for gid in range(2, group_id + 1):
                group_rows = gdf[gdf["tgroup"] == gid]
                tcount = group_rows["dcount"].sum()
                tmean = group_rows["dmean"].iloc[0]
                tvar = group_rows["dvar"].iloc[0]
                dfree = fcount + tcount - 2
                pvar = (
                    ((fcount - 1) * fvar + (tcount - 1) * tvar) / dfree
                    if dfree > 0
                    else fvar
                )
                t_stat = (fmean - tmean) / np.sqrt(pvar * (1 / fcount + 1 / tcount))
                p_val = 1 - student_t.cdf(abs(t_stat), dfree)
                if p_val <= ftalpha:
                    fg += 1
                    fcount, fmean, fvar = tcount, tmean, tvar
                else:
                    new_count = fcount + tcount
                    new_mean = (fcount * fmean + tcount * tmean) / new_count
                    new_var = (
                        ((fcount - 1) * fvar + (tcount - 1) * tvar) / (new_count - 2)
                        if new_count > 2
                        else fvar
                    )
                    fcount, fmean, fvar = new_count, new_mean, new_var
                for _, grp in group_rows.iterrows():
                    final_assign[grp["value"]] = fg

            # 4. Write recode file

            fa = (
                pd.DataFrame.from_dict(final_assign, orient="index", columns=["fgroup"])
                .rename_axis(var)
                .reset_index()
                .merge(tmp[[var, "dcount", "dmean"]], on=var)
            )

            grouped = fa.groupby("fgroup").agg(
                total_count=("dcount", "sum"),
                weighted_sum=("dmean", lambda x: (x * fa.loc[x.index, "dcount"]).sum()),
            )
            grouped["dmean"] = grouped["weighted_sum"] / grouped["total_count"]
            grouped["INDEX"] = grouped["dmean"] / overall_avg * 100

            group_stats = grouped[["dmean", "INDEX"]].reset_index()

            val_col = "INDEX" if nom_method.upper() == "INDEX" else "dmean"
            mapping = {
                row[var]: row[val_col]
                for _, row in fa[["fgroup", var]]
                .merge(group_stats, on="fgroup")
                .iterrows()
            }

            os.makedirs(path_output, exist_ok=True)
            recode_py = os.path.join(path_output, "CE2_Nominal_Var_Recode.py")
            with open(recode_py, "a", encoding="utf-8") as f:
                f.write(f"# Recode variable: {var}\n")
                if nom_method.upper() == "BINARY":
                    for g in range(1, group_stats["fgroup"].max() + 1):
                        vals = [v for v, grp in final_assign.items() if grp == g]
                        f.write(
                            f"df['{prefix}{var}_X{g}'] = "
                            f"df['{var}'].isin({vals}).astype(int)\n"
                        )
                else:
                    f.write(
                        f"df['{prefix}{var}'] = df['{var}'].map({mapping}).fillna(0)\n"
                    )
            if logger:
                logger.info(f"Recode for variable {var} written to {recode_py}")
                logger.debug(f"Recode for variable {var} written to {recode_py}")

            # Update vars2_nom
            cond = vars2_nom["var_name"] == var
            vars2_nom.loc[cond, ["uniq_cnt", "miss_cnt", "max_cnt"]] = [
                uniq_cnt,
                miss_cnt,
                max_cnt,
            ]

            if logger:
                logger.info(
                    f"Updated vars2_nom for {var}: uniq_cnt={uniq_cnt}, miss_cnt={miss_cnt}, max_cnt={max_cnt}"
                )
                logger.debug(
                    f"Updated vars2_nom for {var}: uniq_cnt={uniq_cnt}, miss_cnt={miss_cnt}, max_cnt={max_cnt}"
                )

            # Profiling
            if profiling:
                tmp = df[[dep_var, var]].copy()
                profiles = [
                    {
                        "Variable": var,
                        "Category": "Overall",
                        "Average_DV": overall_avg,
                        "Count": nobs,
                        "Percent": 1.0,
                        "Index": 100.0,
                        "Star": "",
                    }
                ]
                group_vals = {}
                for val, grp in final_assign.items():
                    group_vals.setdefault(grp, []).append(val)
                for grp, vals in group_vals.items():
                    cnt = df[df[var].isin(vals)][dep_var].count()
                    avg = df[df[var].isin(vals)][dep_var].mean()
                    pct = cnt / nobs
                    idx = (avg / overall_avg) * 100
                    star = (
                        "* (+)"
                        if idx >= 110
                        else (
                            "  (+)"
                            if idx > 100
                            else (
                                "  (-)"
                                if idx > 90
                                else "* (-)" if idx <= 90 else "  (0)"
                            )
                        )
                    )
                    cat = ",".join(map(str, vals))
                    profiles.append(
                        {
                            "Variable": var,
                            "Category": cat,
                            "Average_DV": avg,
                            "Count": cnt,
                            "Percent": pct,
                            "Index": idx,
                            "Star": star,
                        }
                    )
                prof_nom = pd.DataFrame(profiles)
                if logger:
                    logger.info(f"Profiling for variable {var} completed.")
                    logger.debug(f"Profiling for variable {var} completed.")
    return vars2_nom, prof_nom


def nom_cntl(
    df: pd.DataFrame,
    nom_vars: list,
    vars2_nom: pd.DataFrame,
    typ_map: dict,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Controller for processing a list of nominal variables.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    nom_vars : list
        Names of nominal variables to process.
    vars2_nom : pd.DataFrame
        Metadata table for nominal variables.
    typ_map : dict
        Mapping of variable names to their data type.
    config : dict
        Configuration dictionary passed to :func:`pnom`.
    logger : logging.Logger, optional
        Optional logger for progress output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Updated ``vars2_nom``, profiling DataFrame and keep list for recoded
        variables.
    """

    path_output = config["path_output"]
    prefix = config.get("prefix", "R1_")
    dep_var = config["dep_var"]

    if logger:
        logger.info(f"nom_cntl_start, number of nominal variables = {len(nom_vars)}")
        logger.debug(f"nom_cntl_start, number of nominal variables = {len(nom_vars)}")
    #
    os.makedirs(path_output, exist_ok=True)
    recode_py = os.path.join(path_output, "CE2_Nominal_Var_Recode.py")
    header = ["# -*- coding: utf-8 -*-", "# Auto-generated nominal recode", ""]
    with open(recode_py, "w") as f:
        f.write("\n".join(header) + "\n")
    #
    prof_nom_all = []
    new_vars = []
    nobs = len(df)
    overall_avg = df[dep_var].mean()
    #
    for var in nom_vars:
        vars2_nom, prof_nom = pnom(
            df=df,
            var=var,
            typ=typ_map[var],
            dep_var=dep_var,
            vars2_nom=vars2_nom,
            overall_avg=overall_avg,
            nobs=nobs,
            config=config,
            logger=logger,
        )
        if prof_nom is not None:
            prof_nom_all.append(prof_nom)
        vars2_nom.loc[vars2_nom["var_name"] == var, "new_var"] = ""
        if vars2_nom.loc[vars2_nom["var_name"] == var, "overall_good"].iloc[0] == 1:
            vars2_nom.loc[vars2_nom["var_name"] == var, "new_var"] = prefix + var
            new_vars.append(prefix + var)
        if logger:
            logger.info(
                f"Processed variable: {var}, overall_good={vars2_nom.loc[vars2_nom['var_name'] == var, 'overall_good'].iloc[0]}"
            )
            logger.debug(
                f"Processed variable: {var}, overall_good={vars2_nom.loc[vars2_nom['var_name'] == var, 'overall_good'].iloc[0]}"
            )
    #
    with open(recode_py, "a") as f:
        f.write("\n# KEEP_LIST_N\n")
        f.write(f"KEEP_LIST_N = {new_vars}\n")

    keep_vars_nom = vars2_nom[vars2_nom["overall_good"] == 1][
        ["var_name", "new_var"]
    ].copy()
    keep_vars_nom = keep_vars_nom.rename(
        columns={"var_name": "orig_var", "new_var": "variable"}
    )

    profile_df = (
        pd.concat(prof_nom_all, ignore_index=True) if prof_nom_all else pd.DataFrame()
    )

    if logger:
        logger.info(f"nom_cntl end, number of new nominal variables = {len(new_vars)}")
        logger.debug(f"nom_cntl end, number of new nominal variables = {len(new_vars)}")

    return vars2_nom, profile_df, keep_vars_nom


def pord(
    df: pd.DataFrame,
    var: str,
    vars2_ord: pd.DataFrame,
    overall_avg: float,
    nobs: int,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process an ordinal variable and generate profiling information.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    var : str
        Ordinal variable to process.
    vars2_ord : pd.DataFrame
        Metadata table for ordinal variables.
    overall_avg : float
        Overall average of the dependent variable.
    nobs : int
        Number of observations.
    config : dict
        Configuration dictionary with options such as ``concrate`` and
        ``num_category``.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Updated ``vars2_ord`` and a profiling DataFrame for ``var``.
    """

    dep_var = config["dep_var"]
    concrate = config.get("concrate", 0.9)
    profiling = config.get("profiling", "Y") == "Y"
    num_category = config.get("num_category", 10)

    if logger:
        logger.info(f"Processing ordinal variable: {var} in pord")
        logger.debug(f"Processing ordinal variable: {var} in pord")
    prof_ord = pd.DataFrame()

    # 1. Summarize counts
    cnts = df[var].value_counts(dropna=False)
    max_cnt = int(cnts.max())
    uniq_cnt = cnts.shape[0]
    # 2. Check constant dep_var in nonmissing
    nonmiss = df.loc[df[var].notna(), dep_var]
    cons_dep = nonmiss.nunique() == 1
    #
    concrate_upper = nobs * concrate
    concrate_good = concrate_upper >= max_cnt
    #
    if not concrate_good:
        if logger:
            logger.warning(
                f"Variable {var} is too concentrated to use: max_cnt {max_cnt} > concrate_upper {concrate_upper}]"
            )
            logger.debug(
                f"Variable {var} is too concentrated to use: max_cnt {max_cnt} > concrate_upper {concrate_upper}]"
            )
        return vars2_ord, prof_ord
    else:
        if cons_dep:
            if logger:
                logger.warning(
                    f"Variable {var} has constant dependent value in nonmissing part and will be excluded"
                )
                logger.debug(
                    f"Variable {var} has constant dependent value in nonmissing part and will be excluded"
                )
            return vars2_ord, prof_ord
        else:
            vars2_ord = pnum(
                df=df, var=var, typ="O", vars2=vars2_ord, config=config, logger=logger
            )
            cond = vars2_ord["var_name"] == var
            vars2_ord.loc[cond, ["overall_good"]] = 1
            if logger:
                logger.info(f"Updated vars2_ord for {var}")
                logger.debug(f"Updated vars2_ord for {var}")
            if profiling:
                if uniq_cnt <= num_category:
                    tmp_prof = prof1(df=df, var=var, config=config, logger=logger)
                else:
                    tmp_prof = prof2(df=df, var=var, config=config, logger=logger)
                prof3_df = prof3(
                    prof=tmp_prof,
                    var=var,
                    overall_avg=overall_avg,
                    nobs=nobs,
                    config=config,
                    logger=logger,
                )
                prof_ord = pd.DataFrame(prof3_df)
            if logger:
                logger.info(f"Profiling for ordinal variable {var} completed.")
                logger.debug(f"Profiling for ordinal variable {var} completed.")
    return vars2_ord, prof_ord


def ord_cntl(
    df: pd.DataFrame,
    ord_vars: list,
    vars2_ord: pd.DataFrame,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Controller for processing a list of ordinal variables.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    ord_vars : list
        Names of ordinal variables to process.
    vars2_ord : pd.DataFrame
        Metadata table for ordinal variables.
    config : dict
        Configuration dictionary passed to :func:`pord`.
    logger : logging.Logger, optional
        Optional logger for progress output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Updated ``vars2_ord``, profiling DataFrame and keep list for recoded
        variables.
    """

    path_output = config["path_output"]
    prefix = config.get("prefix", "R1_")
    dep_var = config["dep_var"]
    missrate = config.get("missrate", 0.75)

    if logger:
        logger.info(f"ord_cntl start, number of ordinal variables = {len(ord_vars)}")
        logger.debug(f"ord_cntl start, number of ordinal variables = {len(ord_vars)}")
    #
    os.makedirs(path_output, exist_ok=True)
    recode_py = os.path.join(path_output, "CE2_Ordinal_Var_Recode.py")
    header = ["# -*- coding: utf-8 -*-", "# Auto-generated ordinal recode", ""]
    with open(recode_py, "w") as f:
        f.write("\n".join(header) + "\n")
    #
    prof_ord_all = []
    new_vars = []
    nobs = len(df)
    overall_avg = df[dep_var].mean()
    #
    for col in ("miss_cnt", "uniq_cnt", "max_cnt", "new_var", "Sign", "Relationship"):
        if col not in vars2_ord.columns:
            vars2_ord[col] = np.nan
    #
    for var in ord_vars:
        miss_cnt = int(df[var].isna().sum())
        uniq_cnt = int(df[var].nunique(dropna=False))
        max_cnt = int(df[var].value_counts(dropna=False).max())
        #
        m = vars2_ord["var_name"] == var
        vars2_ord.loc[m, "miss_cnt"] = miss_cnt
        vars2_ord.loc[m, "uniq_cnt"] = uniq_cnt
        vars2_ord.loc[m, "max_cnt"] = max_cnt
        #
        if miss_cnt > nobs * missrate:
            if logger:
                logger.warning(
                    f"Variable {var} missing count is {miss_cnt}, which is too high"
                )
                logger.debug(
                    f"Variable {var} missing count is {miss_cnt}, which is too high"
                )
        else:
            vars2_ord, prof_ord = pord(
                df=df,
                var=var,
                vars2_ord=vars2_ord,
                overall_avg=overall_avg,
                nobs=nobs,
                config=config,
                logger=logger,
            )
            if prof_ord is not None:
                prof_ord_all.append(prof_ord)
            vars2_ord.loc[vars2_ord["var_name"] == var, "new_var"] = ""
            if vars2_ord.loc[vars2_ord["var_name"] == var, "overall_good"].iloc[0] == 1:
                vars2_ord.loc[vars2_ord["var_name"] == var, "new_var"] = prefix + var
                new_vars.append(prefix + var)
            if logger:
                logger.info(
                    f"Processed variable: {var}, overall_good={vars2_ord.loc[vars2_ord['var_name'] == var, 'overall_good'].iloc[0]}"
                )
                logger.debug(
                    f"Processed variable: {var}, overall_good={vars2_ord.loc[vars2_ord['var_name'] == var, 'overall_good'].iloc[0]}"
                )
    #
    with open(recode_py, "a") as f:
        f.write("\n# KEEP_LIST_O\n")
        f.write(f"KEEP_LIST_O = {new_vars}\n")

    keep_vars_ord = vars2_ord[vars2_ord["overall_good"] == 1][
        ["var_name", "new_var"]
    ].copy()
    keep_vars_ord = keep_vars_ord.rename(
        columns={"var_name": "orig_var", "new_var": "variable"}
    )

    profile_df = (
        pd.concat(prof_ord_all, ignore_index=True) if prof_ord_all else pd.DataFrame()
    )

    if logger:
        logger.info(f"ord_cntl end, number of new ordinal variables = {len(new_vars)}")
        logger.debug(f"ord_cntl end, number of new ordinal variables = {len(new_vars)}")

    return vars2_ord, profile_df, keep_vars_ord


def pcont(
    df: pd.DataFrame,
    var: str,
    vars2_cont: pd.DataFrame,
    overall_avg: float,
    nobs: int,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Process a continuous variable using :func:`pnum` and profiling.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    var : str
        Continuous variable to process.
    vars2_cont : pd.DataFrame
        Metadata table for continuous variables.
    overall_avg : float
        Overall mean of the dependent variable.
    nobs : int
        Number of observations.
    config : dict
        Configuration dictionary.
    logger : logging.Logger, optional
        Optional logger for debug output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        Updated ``vars2_cont`` and profiling DataFrame for ``var``.
    """

    dep_var = config["dep_var"]
    profiling = config.get("profiling", "Y") == "Y"

    if logger:
        logger.info(f"Processing continuous variable: {var} in pcont")
        logger.debug(f"Processing continuous variable: {var} in pcont")
    prof_cont = pd.DataFrame()

    nonmiss = df.loc[df[var].notna(), dep_var]
    cons_dep = nonmiss.nunique() == 1
    if cons_dep:
        if logger:
            logger.warning(
                f"Variable {var} has constant dependent value in nonmissing part and will be excluded"
            )
            logger.debug(
                f"Variable {var} has constant dependent value in nonmissing part and will be excluded"
            )
        return vars2_cont, prof_cont
    else:
        vars2_cont = pnum(
            df=df, var=var, typ="C", vars2=vars2_cont, config=config, logger=logger
        )
        cond = vars2_cont["var_name"] == var
        vars2_cont.loc[cond, ["overall_good"]] = 1
        if logger:
            logger.info(f"Updated vars2_cont for {var}")
            logger.debug(f"Updated vars2_cont for {var}")
        if profiling:
            tmp_prof = prof2(df=df, var=var, config=config, logger=logger)
            prof3_df = prof3(
                prof=tmp_prof,
                var=var,
                overall_avg=overall_avg,
                nobs=nobs,
                config=config,
                logger=logger,
            )
            prof_cont = pd.DataFrame(prof3_df)
        if logger:
            logger.info(f"Profiling for continuous variable {var} completed.")
            logger.debug(f"Profiling for continuous variable {var} completed.")
    return vars2_cont, prof_cont


def cont_cntl(
    df: pd.DataFrame,
    cont_vars: list,
    vars2_cont: pd.DataFrame,
    config: dict,
    logger: Optional[logging.Logger] = None,
) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    """Controller for processing continuous variables.

    Parameters
    ----------
    df : pd.DataFrame
        Input dataset.
    cont_vars : list
        Names of continuous variables to process.
    vars2_cont : pd.DataFrame
        Metadata table for continuous variables.
    config : dict
        Configuration dictionary passed to :func:`pcont`.
    logger : logging.Logger, optional
        Optional logger for progress output.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
        Updated ``vars2_cont``, profiling DataFrame and keep list for recoded
        variables.
    """

    path_output = config["path_output"]
    prefix = config.get("prefix", "R1_")
    dep_var = config["dep_var"]
    missrate = config.get("missrate", 0.75)
    p_lo = config.get("p_lo", 1)
    p_hi = config.get("p_hi", 99)

    if logger:
        logger.info(
            f"cont_cntl start, number of continuous variables = {len(cont_vars)}"
        )
        logger.debug(
            f"cont_cntl start, number of continuous variables = {len(cont_vars)}"
        )
    #
    os.makedirs(path_output, exist_ok=True)
    recode_py = os.path.join(path_output, "CE2_Continuous_Var_Recode.py")
    header = ["# -*- coding: utf-8 -*-", "# Auto-generated continuous recode", ""]
    with open(recode_py, "w") as f:
        f.write("\n".join(header) + "\n")
    #
    prof_cont_all = []
    new_vars = []
    nobs = len(df)
    overall_avg = df[dep_var].mean()
    #
    for col in (
        "miss_cnt",
        "P_lo",
        "P_hi",
        "overall_good",
        "new_var",
        "Sign",
        "Relationship",
    ):
        if col not in vars2_cont.columns:
            vars2_cont[col] = np.nan
    #
    for var in cont_vars:
        miss_cnt = int(df[var].isna().sum())
        lo = df[var].quantile(p_lo / 100)
        hi = df[var].quantile(p_hi / 100)
        #
        m = vars2_cont["var_name"] == var
        vars2_cont.loc[m, "miss_cnt"] = miss_cnt
        vars2_cont.loc[m, "P_lo"] = lo
        vars2_cont.loc[m, "P_hi"] = hi
        #
        if miss_cnt > nobs * missrate:
            if logger:
                logger.warning(
                    f"Variable {var} missing count is {miss_cnt}, which is too high"
                )
                logger.debug(
                    f"Variable {var} missing count is {miss_cnt}, which is too high"
                )
        else:
            if lo == hi:
                if logger:
                    logger.warning(
                        f"Variable {var} has constant value and will be excluded"
                    )
                    logger.warning(
                        f"Variable {var} has constant value and will be excluded"
                    )
            else:
                vars2_cont, prof_cont = pcont(
                    df=df,
                    var=var,
                    vars2_cont=vars2_cont,
                    overall_avg=overall_avg,
                    nobs=nobs,
                    config=config,
                    logger=logger,
                )
                if prof_cont is not None:
                    prof_cont_all.append(prof_cont)
                vars2_cont.loc[vars2_cont["var_name"] == var, "new_var"] = ""
                if (
                    vars2_cont.loc[vars2_cont["var_name"] == var, "overall_good"].iloc[0]
                    == 1
                ):
                    vars2_cont.loc[vars2_cont["var_name"] == var, "new_var"] = prefix + var
                    new_vars.append(prefix + var)
                if logger:
                    logger.info(
                        f"Processed variable: {var}, overall_good={vars2_cont.loc[vars2_cont['var_name'] == var, 'overall_good'].iloc[0]}"
                    )
                    logger.debug(
                        f"Processed variable: {var}, overall_good={vars2_cont.loc[vars2_cont['var_name'] == var, 'overall_good'].iloc[0]}"
                    )
    #
    with open(recode_py, "a") as f:
        f.write("\n# KEEP_LIST_C\n")
        f.write(f"KEEP_LIST_C = {new_vars}\n")

    keep_vars_cont = vars2_cont[vars2_cont["overall_good"] == 1][
        ["var_name", "new_var"]
    ].copy()
    keep_vars_cont = keep_vars_cont.rename(
        columns={"var_name": "orig_var", "new_var": "variable"}
    )

    profile_df = (
        pd.concat(prof_cont_all, ignore_index=True) if prof_cont_all else pd.DataFrame()
    )

    if logger:
        logger.info(
            f"cont_cntl end, number of new continuous variables = {len(new_vars)}"
        )
        logger.debug(
            f"cont_cntl end, number of new continuous variables = {len(new_vars)}"
        )

    return vars2_cont, profile_df, keep_vars_cont


def CE_EDA_Recode(
    indsn: pd.DataFrame, config: dict, logger: Optional[logging.Logger] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Run the full EDA and recoding workflow.

    Parameters
    ----------
    indsn : pd.DataFrame
        Input dataset.
    config : dict
        Workflow configuration.
    logger : logging.Logger, optional
        Optional logger for progress messages.

    Returns
    -------
    Tuple[pd.DataFrame, pd.DataFrame]
        The recoded dataset and profiling DataFrame.
    """

    if logger:
        logger.info("*** 2. EDA, profiling and Recode ***")
        logger.debug(
            f"Input dataset rows: {len(indsn)} columns: {indsn.columns.tolist()}"
        )

    # Variable Lists
    bin_vars = config.get("bin_vars", [])
    nom_vars = config.get("nom_vars", [])
    ord_vars = config.get("ord_vars", [])
    cont_vars = config.get("cont_vars", [])
    # General Macro Variables
    path_output = config["path_output"]
    id = config["id"]
    dep_var = config["dep_var"]
    # Optional Macro Variables: These all have defaults that can be used
    # General Macro Variables
    # Macro 2: Recoding macro variables
    profiling = config.get("profiling", "Y").upper() == "Y"
    minbinnc = config.get("minbinnc", 500)
    minbinnp = config.get("minbinnp", 0.05)

    # 1. Build workfile (exclude mod_val_test==3)
    mask_val = indsn.get("mod_val_test", 0) == 3
    n_val = mask_val.sum()
    workfile = indsn.loc[mask_val == False].copy()

    if logger:
        logger.info(f"Excluded Validation rows: {n_val}")
        logger.debug(f"Rows after exclusion: {len(workfile)}")

    # 2. Global stats
    nobs = len(workfile)
    overall_avg = workfile[dep_var].mean()
    minbinn = max(minbinnc, int(nobs * minbinnp))

    # 3. Prepare profile DataFrame and var lookup
    profile_df = pd.DataFrame()
    var_lookup = []
    # UAT workfile['IsChild'] = workfile['IsChild'].map({1: 'y', 0: 'N'})
    # UAT workfile.loc[9:29, 'IsChild'] = np.nan
    # UAT workfile['Embarked'] =
    # UAT workfile['Cabin2'] =
    # UAT workfile.loc[9:29, 'Pclass'] = np.nan

    # 4. Binary variables
    if bin_vars:
        vars2_bin = pd.DataFrame({"var_name": bin_vars})
        typ_map = {
            v: ("numeric" if pd.api.types.is_numeric_dtype(workfile[v]) else "str")
            for v in bin_vars
        }
        vars2_bin, prof_bin, keep_vars_bin = bin_cntl(
            df=workfile,
            bin_vars=bin_vars,
            vars2_bin=vars2_bin,
            typ_map=typ_map,
            config=config,
            logger=logger,
        )
        if profiling and not prof_bin.empty:
            profile_df = pd.concat([profile_df, prof_bin], ignore_index=True)
        if logger:
            logger.info(f"Binary variables processed: {len(vars2_bin)}")
            logger.debug(f"Binary variables processed: {len(vars2_bin)}")
        if not keep_vars_bin.empty:
            keep_vars_bin_list = keep_vars_bin["variable"].dropna().tolist()
            vars2_bin["type"] = "binary"
            col = "type"
            cols = list(vars2_bin.columns)
            cols.pop(cols.index(col))
            cols.insert(1, col)
            vars2_bin = vars2_bin[cols]
        else:
            keep_vars_bin_list = []
    else:
        typ_map = {}
        vars2_bin = pd.DataFrame()
        keep_vars_bin_list = []
        if logger:
            logger.info("No binary variables to process.")
            logger.debug("No binary variables to process.")

    # 5. Nominal variables
    if nom_vars:
        vars2_nom = pd.DataFrame({"var_name": nom_vars})
        typ_map = {
            v: ("numeric" if pd.api.types.is_numeric_dtype(workfile[v]) else "str")
            for v in nom_vars
        }
        vars2_nom, prof_nom, keep_vars_nom = nom_cntl(
            df=workfile,
            nom_vars=nom_vars,
            vars2_nom=vars2_nom,
            typ_map=typ_map,
            config=config,
            logger=logger,
        )
        if profiling and not prof_nom.empty:
            profile_df = pd.concat([profile_df, prof_nom], ignore_index=True)
        if logger:
            logger.info(f"Nominal variables processed: {len(vars2_nom)}")
            logger.debug(f"Nominal variables processed: {len(vars2_nom)}")
        if not keep_vars_nom.empty:
            keep_vars_nom_list = keep_vars_nom["variable"].dropna().tolist()
            vars2_nom["type"] = "nominal"
            col = "type"
            cols = list(vars2_nom.columns)
            cols.pop(cols.index(col))
            cols.insert(1, col)
            vars2_nom = vars2_nom[cols]
        else:
            keep_vars_nom_list = []
    else:
        typ_map = {}
        vars2_nom = pd.DataFrame()
        keep_vars_nom_list = []
        if logger:
            logger.info("No nominal variables to process.")
            logger.debug("No nominal variables to process.")

    # 6. Ordinal variables
    if ord_vars:
        vars2_ord = pd.DataFrame({"var_name": ord_vars})
        vars2_ord, prof_ord, keep_vars_ord = ord_cntl(
            df=workfile,
            ord_vars=ord_vars,
            vars2_ord=vars2_ord,
            config=config,
            logger=logger,
        )
        if profiling and not prof_ord.empty:
            profile_df = pd.concat([profile_df, prof_ord], ignore_index=True)
        if logger:
            logger.info(f"Ordinal variables processed: {len(vars2_ord)}")
            logger.debug(f"Ordinal variables processed: {len(vars2_ord)}")
        if not keep_vars_ord.empty:
            keep_vars_ord_list = keep_vars_ord["variable"].dropna().tolist()
            vars2_ord["type"] = "ordinal"
            col = "type"
            cols = list(vars2_ord.columns)
            cols.pop(cols.index(col))
            cols.insert(1, col)
            vars2_ord = vars2_ord[cols]
        else:
            keep_vars_ord_list = []
    else:
        vars2_ord = pd.DataFrame()
        keep_vars_ord_list = []
        if logger:
            logger.info("No Ordinal variables to process.")
            logger.debug("No Ordinal variables to process.")

    # 7. Continuous variables
    if cont_vars:
        vars2_cont = pd.DataFrame({"var_name": cont_vars})
        vars2_cont, prof_cont, keep_vars_cont = cont_cntl(
            df=workfile,
            cont_vars=cont_vars,
            vars2_cont=vars2_cont,
            config=config,
            logger=logger,
        )
        if profiling and not prof_cont.empty:
            profile_df = pd.concat([profile_df, prof_cont], ignore_index=True)
        if logger:
            logger.info(f"Continuous variables processed: {len(vars2_cont)}")
            logger.debug(f"Continuous variables processed: {len(vars2_cont)}")
        if not vars2_cont.empty:
            keep_vars_cont_list = keep_vars_cont["variable"].dropna().tolist()
            vars2_cont["type"] = "continuous"
            col = "type"
            cols = list(vars2_cont.columns)
            cols.pop(cols.index(col))
            cols.insert(1, col)
            vars2_cont = vars2_cont[cols]
        else:
            keep_vars_cont_list = []
    else:
        vars2_cont = pd.DataFrame()
        keep_vars_cont_list = []
        if logger:
            logger.info("No Continuous variables to process.")
            logger.debug("No Continuous variables to process.")

    # 8. Apply recode to entire dataset (ensuring new_vars exist)
    recoded = workfile.copy()
    for script in glob.glob(os.path.join(path_output, "CE2_*_Var_Recode.py")):
        # Load the script
        code = open(script, "r", encoding="utf-8").read()
        # Execute on the recoded DataFrame
        # Ensure the DataFrame name matches and np is available
        local_vars = {"df": recoded, "np": np, "nan": np.nan}
        exec(code, local_vars)
        # Update recoded so new variables are included
        recoded = local_vars["df"]
        exec(code, {"df": recoded, "np": np, "nan": np.nan})

    # Now all new variables can be selected
    keep_list = config.get("keep_list", [])
    base_keep = [id, "mod_val_test", dep_var]
    all_new = (
        keep_vars_bin_list
        + keep_vars_nom_list
        + keep_vars_ord_list
        + keep_vars_cont_list
    )
    CE2_Recoded = recoded[base_keep + keep_list + all_new].copy()

    # 4. Generate summaries for each variable type and output to separate Excel sheets
    # report_path = os.path.join(path_output, "CE2_EDA_report.xlsx")
    # with pd.ExcelWriter(report_path, engine="xlsxwriter") as writer:
    #     # 4.1 Binary
    #     if not vars2_bin.empty:
    #         vars2_bin[vars2_bin["overall_good"] == 1].to_excel(
    #             writer, sheet_name="Binary Variables", index=False
    #         )
    #     # 4.2 Nominal
    #     if not vars2_nom.empty:
    #         vars2_nom[vars2_nom["overall_good"] == 1].to_excel(
    #             writer, sheet_name="Nominal Variables", index=False
    #         )
    #     # 4.3 Ordinal
    #     if not vars2_ord.empty:
    #         vars2_ord[vars2_ord["overall_good"] == 1].to_excel(
    #             writer, sheet_name="Ordinal Variables", index=False
    #         )
    #     # 4.4 Continuous
    #     if not vars2_cont.empty:
    #         vars2_cont[vars2_cont["overall_good"] == 1].to_excel(
    #             writer, sheet_name="Continuous Variables", index=False
    #         )
    #     # 4.5 Correlation
    #     if dep_var in CE2_Recoded.columns and pd.api.types.is_numeric_dtype(
    #         CE2_Recoded[dep_var]
    #     ):
    #         corr = (
    #             CE2_Recoded.corr()[dep_var]
    #             .drop(dep_var)
    #             .sort_values(key=lambda s: s.abs(), ascending=False)
    #         )
    #         corr_df = corr.reset_index().rename(
    #             columns={"index": "variable", dep_var: "correlation"}
    #         )
    #     else:
    #         corr_df = pd.DataFrame(columns=["variable", "correlation"])
    #     corr_df.to_excel(writer, sheet_name="Correlations", index=False)

    return CE2_Recoded, profile_df
