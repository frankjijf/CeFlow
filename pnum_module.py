import pandas as pd
import numpy as np
import statsmodels.api as sm
from typing import Dict


def pnum(
    df: pd.DataFrame,
    var: str,
    dep_var: str,
    binary_dv: bool = True,
    impmethod: str = 'median',
    transformation: bool = True,
    cap_floor: bool = True,
    stdmethod: str = 'NO',
    prefix: str = 'R1_',
    pvalue: float = 0.05,
    min_size: int = 500,
) -> Dict[str, float]:
    """Approximate implementation of the SAS %macro pnum.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data.
    var : str
        Independent variable to process.
    dep_var : str
        Dependent variable name.
    binary_dv : bool, optional
        Whether the dependent variable is binary.
    impmethod : str, optional
        Missing value imputation method, e.g. ``"median"`` or ``"ER"``.
    transformation : bool, optional
        Add SQ_, SR_, LN_ transformed features when True.
    cap_floor : bool, optional
        Apply capping/flooring before transformations.
    stdmethod : str, optional
        Standardization method, ``"STD"`` or ``"NO"``.
    prefix : str, optional
        Prefix for the generated variable name.
    pvalue : float, optional
        P-value threshold for significance (informational).
    min_size : int, optional
        Minimum missing count to apply equal response (ER) method.

    Returns
    -------
    dict
        Dictionary with key statistics and modeling results.
    """
    ser = df[var].dropna()
    if ser.empty:
        raise ValueError(f"{var} has no non-missing values")

    var_mean = ser.mean()
    var_median = ser.median()
    var_min = ser.min()
    var_max = ser.max()
    var_p1, var_p25, var_p75, var_p99 = np.percentile(ser, [1, 25, 75, 99])

    iqr = max(var_p75 - var_p25, var_p99 - var_p75, var_p25 - var_p1)
    var_lb = min(max(var_p25 - 1.5 * iqr, var_min), var_p1)
    var_ub = max(min(var_p75 + 1.5 * iqr, var_max), var_p99)
    if var_lb == var_ub:
        var_lb, var_ub = var_min, var_max
    var_mid = (var_max + var_min) / 2

    method = impmethod.upper()
    if method in {"MEAN", "STD"}:
        var_miss = var_mean
    elif method in {"MEDIAN", "IQR", "MAD"}:
        var_miss = var_median
    elif method == "RANGE":
        var_miss = var_min
    elif method == "MIDRANGE":
        var_miss = var_mid
    elif method in {"SUM", "EUCLEN", "USTD", "MAXABS"}:
        var_miss = 0.0
    else:
        var_miss = var_median

    df_mod = df.copy()
    if cap_floor:
        df_mod[var] = df_mod[var].clip(var_lb, var_ub)

    mod_vars = [var]
    if transformation:
        df_mod[f"SQ_{var}"] = df_mod[var] ** 2
        df_mod[f"SR_{var}"] = np.sqrt(np.maximum(df_mod[var], 0))
        df_mod[f"LN_{var}"] = np.log(np.maximum(df_mod[var], 1e-5))
        mod_vars.extend([f"SQ_{var}", f"SR_{var}", f"LN_{var}"])

    reg_df = df_mod.dropna(subset=[var, dep_var])
    X = sm.add_constant(reg_df[mod_vars])
    y = reg_df[dep_var]
    if binary_dv:
        model = sm.Logit(y, X).fit(disp=False)
    else:
        model = sm.OLS(y, X).fit()

    estimate = model.params.get(var, np.nan)
    intercept = model.params.get('const', np.nan)
    prob = model.pvalues.get(var, np.nan)
    sign = '(+)' if estimate >= 0 else '(-)'
    relationship = sign

    if method == 'ER' and df[var].isna().sum() >= min_size and not np.isnan(estimate) and estimate != 0:
        miss_rr = df.loc[df[var].isna(), dep_var].mean()
        if binary_dv:
            miss_impute = (np.log(miss_rr / (1 - miss_rr)) - intercept) / estimate
        else:
            miss_impute = (miss_rr - intercept) / estimate
        miss_impute = np.clip(miss_impute, var_lb, var_ub)
        var_miss = miss_impute

    if stdmethod.upper() == 'STD':
        loc = reg_df[var].mean()
        scale = reg_df[var].std(ddof=0)
    else:
        loc = None
        scale = None

    new_var = prefix + var

    return {
        'var_LB': var_lb,
        'var_UB': var_ub,
        'var_median': var_median,
        'var_miss': var_miss,
        'new_var': new_var,
        'sign': sign,
        'relationship': relationship,
        'pvalue': prob,
        'intercept': intercept,
        'estimate': estimate,
        'location': loc,
        'scale': scale,
    }


if __name__ == "__main__":
    df = pd.read_csv("Titanic-Dataset.csv")
    result = pnum(df, var="Age", dep_var="Survived", binary_dv=True)
    print(result)
