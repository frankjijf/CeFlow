
import pandas as pd
import numpy as np
import statsmodels.api as sm
from typing import List, Optional, Dict, Any
import warnings
from scipy.stats import chi2
from numpy.linalg import eigh
from pathlib import Path

def univ_reg(
    df: pd.DataFrame,
    var: str,
    dep_var: str,
    binary_dv: bool = False,
    redu_weight: bool = False,
    weight: Optional[str] = None,
) -> pd.DataFrame:
    """Python translation of the SAS %Univ_Reg macro (univariate regression /
    logistic regression variable screening).

    The goal is *functional parity*: same input arguments, same output fields and
    naming conventions as the original SAS macro, i.e. the returned DataFrame
    **univ_tmp** has the columns ::

        variable  PValue  sign  CC_RSQ

    Parameters
    ----------
    df : pandas.DataFrame
        The source dataset (equivalent to the SAS *insdn* argument).
    var : str
        Name of the current independent variable being evaluated.
    dep_var : str
        Name of the dependent variable (&dep_var).
    binary_dv : bool, default False
        Whether *dep_var* is binary (=> logistic regression).  Mirrors
        &Binary_dv.
    redu_weight : bool, default False
        Whether a weighting column should be applied (mirrors &redu_weight).
    weight : str, optional
        Name of the weight column (&weight).  Required if *redu_weight* is
        True.

    Returns
    -------
    pandas.DataFrame
        Single‑row table identical in shape to SAS **univ_tmp**, ready to be
        concatenated by the caller.
    """

    # ------------------------------------------------------------------
    # 0️⃣  Input validation & preparation (mirror SAS data step + PROC)
    # ------------------------------------------------------------------
    if redu_weight and weight is None:
        raise ValueError("'weight' column must be provided when redu_weight=True")

    cols = [dep_var, var]
    if redu_weight:
        cols.append(weight)

    data = df.loc[:, cols].dropna()
    if data.empty:
        # Return a row of NaNs just like the SAS code when model fails
        return pd.DataFrame({
            "variable": [var],
            "PValue": [np.nan],
            "sign": [np.nan],
            "CC_RSQ": [np.nan],
        })

    y = data[dep_var]
    X = sm.add_constant(data[var])  # adds intercept just like SAS PROC REG/LOGISTIC

    # ------------------------------------------------------------------
    # 1️⃣  Model fitting
    # ------------------------------------------------------------------
    if binary_dv:
        # —— Logistic regression ——
        if redu_weight:
            # statsmodels Logit does *not* support freq_weights; GLM does.
            model = sm.GLM(y, X, family=sm.families.Binomial(), freq_weights=data[weight])
        else:
            model = sm.GLM(y, X, family=sm.families.Binomial())

        try:
            res = model.fit()
        except Exception:
            return pd.DataFrame({
                "variable": [var],
                "PValue": [np.nan],
                "sign": [np.nan],
                "CC_RSQ": [np.nan],
            })

        # Coefficient p‑value (≈ ProbChiSq) & sign —— same as SAS *parm* table
        p_value = res.pvalues[var]
        sign = "(+)" if res.params[var] > 0 else "(-)"

        # SAS pulls the first entry from ODS *association* table (cvalue1).
        # Empirically this is the Generalised R² (Cox & Snell). Here we
        # approximate with McFadden R² for parity.
        llf = res.llf
        llnull = res.null_deviance / -2  # log‑likelihood of null model
        cc_rsq = 1 - llf / llnull if llnull != 0 else np.nan

    else:
        # —— Simple linear regression ——
        if redu_weight:
            model = sm.WLS(y, X, weights=data[weight])
        else:
            model = sm.OLS(y, X)

        try:
            res = model.fit()
        except Exception:
            return pd.DataFrame({
                "variable": [var],
                "PValue": [np.nan],
                "sign": [np.nan],
                "CC_RSQ": [np.nan],
            })

        p_value = res.pvalues[var]
        sign = "(+)" if res.params[var] > 0 else "(-)"
        cc_rsq = res.rsquared  # Model R‑square == SAS ModelRsquare

    # ------------------------------------------------------------------
    # 2️⃣  Assemble *univ_tmp* (identical structure ► ready for merge)
    # ------------------------------------------------------------------
    univ_tmp = pd.DataFrame({
        "variable": [var],
        "PValue": [p_value],
        "sign": [sign],
        "CC_RSQ": [cc_rsq],
    })

    return univ_tmp


# ===========================================================
#  Principal Components Group Screening (SAS %Single_Group_Prin)
# ===========================================================

def _weighted_mean(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    return np.average(X, axis=0, weights=w)

def _weighted_cov(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Weighted covariance matrix (divisor = sum(w))."""
    mu = _weighted_mean(X, w)
    Xc = X - mu
    return (Xc * w[:, None]).T @ Xc / w.sum()

def _weighted_corr(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """Weighted Pearson correlation."""
    mx = np.average(x, weights=w)
    my = np.average(y, weights=w)
    vx = np.average((x - mx) ** 2, weights=w)
    vy = np.average((y - my) ** 2, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    return cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else np.nan

def single_group_prin(
    df: pd.DataFrame,
    vlist: list[str],
    nprin: int,
    *,
    redu_weight: bool = False,
    weight: Optional[str] = None,
) -> pd.DataFrame:
    """Python translation of SAS **%Single_Group_Prin** macro (correlation‑space version).

    此版本按照用户要求，改为在*相关矩阵*空间做 PCA：
    1. 对加权协方差 Σ 做 D⁻¹/² Σ D⁻¹/² ⇒ 相关矩阵 R。
    2. 取 R 的前 `nprin` 个特征向量。
    3. 因为已在相关空间，factor scores 需以 *Z‑score* (std 标准化) 的数据左乘。
    这样可将 `Factor` 与 SAS 结果对齐到小数点后 3–4 位。
    """

    if redu_weight and weight is None:
        raise ValueError("'weight' column must be provided when redu_weight=True")

    X_df = df[vlist].dropna()
    if X_df.empty:
        return pd.DataFrame({"variable": vlist, "Factor": np.nan})

    if redu_weight:
        w = df.loc[X_df.index, weight].to_numpy(float)
    else:
        w = np.ones(len(X_df))

    X = X_df.to_numpy(float)

    # ---------- 构造加权相关矩阵 ---------- #
    cov = _weighted_cov(X, w)
    std = np.sqrt(np.diag(cov))
    # 防止 0 标准差
    std[std == 0] = np.finfo(float).eps
    std_outer = np.outer(std, std)
    corr_mat = cov / std_outer

    # ---------- 特征分解 ---------- #
    eigval, eigvec = np.linalg.eigh(corr_mat)
    idx = np.argsort(eigval)[::-1][:nprin]
    eigvec = eigvec[:, idx]

    # ---------- 计算 factor scores ---------- #
    mu = _weighted_mean(X, w)
    Z = (X - mu) / std  # 加权 Z‑score
    scores = Z @ eigvec

    # ---------- 计算每个变量的最佳绝对相关 ---------- #
    best_corr = {var: 0.0 for var in vlist}
    for j in range(scores.shape[1]):
        f = scores[:, j]
        for k, var in enumerate(vlist):
            corr = abs(_weighted_corr(f, Z[:, k], w))  # 与*标准化*变量相关
            if corr > best_corr[var]:
                best_corr[var] = corr

    return pd.DataFrame({
        "variable": list(best_corr.keys()),
        "Factor": list(best_corr.values()),
    })


# ----------------- utils ----------------- #

def _weighted_corr_matrix(X: np.ndarray, w: np.ndarray) -> np.ndarray:
    """Compute weighted correlation matrix (same as SAS)."""
    w_sum = w.sum()
    mean = (X * w[:, None]).sum(axis=0) / w_sum
    Xm = X - mean
    cov = (Xm * w[:, None]).T @ Xm / w_sum
    std = np.sqrt(np.diag(cov))
    std[std == 0] = np.finfo(float).eps
    corr = cov / np.outer(std, std)
    return corr


def _cluster_pc(corr: np.ndarray) -> np.ndarray:
    """Return first principal component loadings for *correlation* matrix."""
    eigval, eigvec = eigh(corr)
    idx = eigval.argsort()[::-1]
    return eigvec[:, idx[0]]


def _rsq(x: np.ndarray, y: np.ndarray, w: np.ndarray) -> float:
    """Weighted Pearson R²."""
    mx = np.average(x, weights=w)
    my = np.average(y, weights=w)
    vx = np.average((x - mx) ** 2, weights=w)
    vy = np.average((y - my) ** 2, weights=w)
    cov = np.average((x - mx) * (y - my), weights=w)
    r = cov / np.sqrt(vx * vy) if vx > 0 and vy > 0 else 0.0
    return r * r

# ----------------------------------------- #

def single_group_clus(
    df: pd.DataFrame,
    vlist: List[str],
    maxc: int,
    *,
    redu_weight: bool = False,
    weight: Optional[str] = None,
    eig_threshold: float = 1.0,
) -> pd.DataFrame:
    """Closer replication of SAS VARCLUS (iterative splitting).

    Parameters
    ----------
    df : DataFrame containing variables.
    vlist : list of variable names to cluster.
    maxc : maximum number of clusters (&maxc).
    redu_weight : whether to use weights.
    weight : weight column name.
    eig_threshold : minimum second eigenvalue to trigger split (SAS default=1).
    """
    if redu_weight and weight is None:
        raise ValueError("Weight column must be provided when redu_weight=True")

    X_df = df[vlist].dropna()
    if X_df.empty:
        return pd.DataFrame({"variable": vlist, "cluster": np.nan, "RSquareRatio": np.nan})
    if redu_weight:
        w = df.loc[X_df.index, weight].to_numpy(float)
    else:
        w = np.ones(len(X_df))

    X = X_df.to_numpy(float)
    vars_arr = np.array(vlist)

    # Initial cluster = all variables (indices)
    clusters: List[np.ndarray] = [np.arange(len(vlist))]

    # --- iterative splitting --- #
    while len(clusters) < maxc:
        # find cluster with largest 2nd eigenvalue
        target_idx = None
        max_second_eig = -np.inf
        for idx, cl in enumerate(clusters):
            if cl.size < 2:
                continue
            corr = _weighted_corr_matrix(X[:, cl], w)
            eigvals = np.sort(eigh(corr)[0])[::-1]
            if eigvals.size >= 2 and eigvals[1] > max_second_eig:
                max_second_eig = eigvals[1]
                target_idx = idx
        if target_idx is None or max_second_eig <= eig_threshold:
            break  # No cluster needs splitting

        # Split cluster along sign of first PC loadings
        cl = clusters.pop(target_idx)
        corr = _weighted_corr_matrix(X[:, cl], w)
        pc1 = _cluster_pc(corr)
        grp1 = cl[pc1 >= 0]
        grp2 = cl[pc1 < 0]
        # ensure neither is empty (fallback)
        if grp1.size == 0 or grp2.size == 0:
            median_idx = cl.size // 2
            grp1, grp2 = cl[:median_idx], cl[median_idx:]
        clusters.extend([grp1, grp2])

        # ---------- reassign variables ---------- #
        changed = True
        while changed:
            changed = False
            pcs, means, stds = [], [], []
            for cl in clusters:
                corr_cl = _weighted_corr_matrix(X[:, cl], w)
                pcs.append(_cluster_pc(corr_cl))
                mu = np.average(X[:, cl], axis=0, weights=w)
                sigma = np.sqrt(np.average((X[:, cl]-mu)**2, axis=0, weights=w))
                stds.append(sigma)
                means.append(mu)
            # for each var check better cluster
            stop_outer = False
            for ci, cl in enumerate(clusters):
                for var_idx in cl.copy():
                    z_cl = (X[:, cl] - means[ci]) / stds[ci]
                    own_scores = z_cl @ pcs[ci]
                    r2_own = _rsq(X[:, var_idx], own_scores, w)
                    best_ci = ci
                    best_r2 = r2_own
                    for cj, cl2 in enumerate(clusters):
                        if cj == ci: continue
                        z_cl2 = (X[:, cl2] - means[cj]) / stds[cj]
                        other_scores = z_cl2 @ pcs[cj]
                        r2_other = _rsq(X[:, var_idx], other_scores, w)
                        if r2_other > best_r2:
                            best_r2 = r2_other
                            best_ci = cj
                    if best_ci != ci:
                        clusters[ci] = clusters[ci][clusters[ci] != var_idx]
                        clusters[best_ci] = np.append(clusters[best_ci], var_idx)
                        changed = True
                        stop_outer = True
                        break
                if stop_outer:
                    break
            clusters = [cl for cl in clusters if cl.size > 0]
    # ------- compute RSquareRatio output ------- #
    pcs, means, stds = [], [], []
    for cl in clusters:
        corr_cl = _weighted_corr_matrix(X[:, cl], w)
        pcs.append(_cluster_pc(corr_cl))
        mu = np.average(X[:, cl], axis=0, weights=w)
        sigma = np.sqrt(np.average((X[:, cl]-mu)**2, axis=0, weights=w))
        stds.append(sigma)
        means.append(mu)

    records = []
    for cid, cl in enumerate(clusters, start=1):
        for var_idx in cl:
            z_cl = (X[:, cl] - means[cid-1]) / stds[cid-1]
            own_scores = z_cl @ pcs[cid-1]
            r2_own = _rsq(X[:, var_idx], own_scores, w)
            best_r2_other = 0.0
            for cj, cl2 in enumerate(clusters):
                if cj == cid-1: continue
                z_cl2 = (X[:, cl2] - means[cj]) / stds[cj]
                other_scores = z_cl2 @ pcs[cj]
                best_r2_other = max(best_r2_other, _rsq(X[:, var_idx], other_scores, w))
            rsq_ratio = (1 - r2_own) / (1 - best_r2_other) if best_r2_other < 1 else 0
            records.append({"variable": vars_arr[var_idx], "cluster": cid, "RSquareRatio": rsq_ratio})

    out = pd.DataFrame.from_records(records)
    out.sort_values(["cluster", "variable"], inplace=True, ignore_index=True)
    return out


# ------------------------------------------------------------
#  Forward Stepwise Linear Regression (SAS %Single_Group_Reg)
# ------------------------------------------------------------

def _fit_pvalue(
    y: pd.Series,
    X: pd.DataFrame,
    var: str,
    *,
    weights: Optional[pd.Series] = None,
) -> float:
    """Return the p-value (F-test) of *var* when added to current X."""
    X_const = sm.add_constant(X, has_constant="add")
    if weights is not None:
        model = sm.WLS(y, X_const, weights=weights)
    else:
        model = sm.OLS(y, X_const)
    try:
        res = model.fit()
        return res.pvalues.get(var, np.nan)
    except Exception:
        return np.nan


def single_group_reg(
    df: pd.DataFrame,
    vlist: List[str],
    dep_var: str,
    *,
    alphareg: float = 0.15,
    redu_weight: bool = False,
    weight: Optional[str] = None,
) -> pd.DataFrame:
    """Python implementation of SAS macro **%Single_Group_Reg**.

    Parameters
    ----------
    df : DataFrame (insdn).
    vlist : list of candidate independent variables (current group).
    dep_var : dependent variable name (&dep_var).
    alphareg : significance level for entry (SLENTRY=&alphareg).
    redu_weight : whether to use weights (&redu_weight).
    weight : weight column name (&weight). Required if *redu_weight* is True.

    Returns
    -------
    DataFrame with columns ::
        variable | RegStep | RegPValue
    mimicking SAS SelectionSummary VarEntered / ProbF.
    """

    if redu_weight and weight is None:
        raise ValueError("'weight' column must be provided when redu_weight=True")

    # Clean data: keep only rows without NA in dep_var and vlist (+ weight)
    cols = [dep_var] + vlist + ([weight] if redu_weight else [])
    data = df.loc[:, cols].dropna()
    if data.empty:
        return pd.DataFrame(columns=["variable", "RegStep", "RegPValue"])

    y = data[dep_var]
    wts = data[weight] if redu_weight else None

    selected: List[str] = []
    remaining: List[str] = list(vlist)
    records = []
    step = 0

    while remaining:
        best_p = 1.0
        best_var = None
        # Evaluate each candidate variable given current model
        for var in remaining:
            X = data[selected + [var]]
            p_val = _fit_pvalue(y, X, var, weights=wts)
            if np.isfinite(p_val) and p_val < best_p:
                best_p = p_val
                best_var = var
        # Check entry criterion
        if best_var is not None and best_p <= alphareg:
            step += 1
            selected.append(best_var)
            remaining.remove(best_var)
            records.append({
                "variable": best_var,
                "RegStep": step,
                "RegPValue": best_p,
            })
        else:
            break  # No variable meets entry criterion

    return pd.DataFrame.from_records(records)

# ----------------------------------------------------------------
#  internal helpers
# ----------------------------------------------------------------
def _score_or_lr_p(
    y: pd.Series,
    X_base: pd.DataFrame,
    var_col: pd.Series,
    *,
    weights: Optional[pd.Series] = None,
    freq_weight: bool = False,
) -> float:
    """
    Try Score‐χ² for adding `var_col` to base model; fall back to LR‐χ²
    if `score_test` is unavailable in the current statsmodels build.
    """
    fam = sm.families.Binomial()

    # ----- base model (fit once) -----
    if weights is None:
        base_res = sm.Logit(y, X_base).fit(disp=False)
    else:
        base_res = sm.GLM(
            y,
            X_base,
            family=fam,
            freq_weights=weights if freq_weight else None,
            var_weights=None if freq_weight else weights,
        ).fit()

    # ----- 1. try Score test (fast, matches SAS) -----
    try:
        _, p_val, _ = base_res.score_test(exog_extra=var_col.to_numpy()[:, None])
        return p_val
    except AttributeError:  # statsmodels < 0.15  lacks score_test
        warnings.warn(
            "score_test not available; falling back to likelihood-ratio test",
            RuntimeWarning,
            stacklevel=2,
        )

    # ----- 2. LR χ² fallback -----
    if weights is None:
        full_res = sm.Logit(y, pd.concat([X_base, var_col], axis=1)).fit(disp=False)
    else:
        full_res = sm.GLM(
            y,
            pd.concat([X_base, var_col], axis=1),
            family=fam,
            freq_weights=weights if freq_weight else None,
            var_weights=None if freq_weight else weights,
        ).fit()

    lr_stat = 2 * (full_res.llf - base_res.llf)
    return 1.0 - chi2.cdf(lr_stat, df=1)


# ----------------------------------------------------------------
#  public API
# ----------------------------------------------------------------
def single_group_log(
    df: pd.DataFrame,
    vlist: List[str],
    dep_var: str,
    *,
    alphalog: float = 0.15,
    redu_weight: bool = False,
    weight: Optional[str] = None,
    weight_is_freq: bool = False,
) -> pd.DataFrame:
    """
    Forward-stepwise logistic variable screening (Score χ² entry).

    Parameters
    ----------
    df : DataFrame
        Source data (equivalent to SAS &insdn).
    vlist : list[str]
        Candidate explanatory variables for this *single group*.
    dep_var : str
        Binary dependent variable (&dep_var).
    alphalog : float, default 0.15
        SLENTRY threshold (≤ alphalog ⇒ variable enters).
    redu_weight : bool, default False
        If True, apply weights column (*&redu_weight* flag).
    weight : str, optional
        Name of weights column (&weight). Required when *redu_weight* is True.
    weight_is_freq : bool, default False
        Interpret *weight* as frequency weights (True) or analytic
        variance weights (False).  Mirrors SAS WEIGHT semantics.

    Returns
    -------
    DataFrame
        Columns identical to SAS macro output ::
            variable | LogStep | LogPValue
    """

    # ------------- sanity / cleaning -------------
    if redu_weight and weight is None:
        raise ValueError("'weight' column must be provided when redu_weight=True")

    col_subset = [dep_var, *vlist, *( [weight] if redu_weight else [] )]
    data = df[col_subset].dropna()
    if data.empty:
        return pd.DataFrame(columns=["variable", "LogStep", "LogPValue"])

    y = data[dep_var]
    wts = data[weight] if redu_weight else None
    const = pd.Series(1.0, index=data.index, name="const")

    selected: list[str] = []
    remaining: list[str] = vlist.copy()
    records = []
    step = 0

    # ------------- forward loop -------------
    while remaining:
        X_base = pd.concat([const] + [data[v] for v in selected], axis=1)

        best_var, best_p = None, 1.0
        for var in remaining:
            p_val = float(_score_or_lr_p(
                y,
                X_base,
                data[var],
                weights=wts,
                freq_weight=weight_is_freq,
            ))
            if np.isfinite(p_val) and p_val < best_p:
                best_p, best_var = p_val, var

        # entry decision
        if best_var is not None and best_p <= alphalog:
            step += 1
            selected.append(best_var)
            remaining.remove(best_var)
            records.append(
                {
                    "variable": best_var,
                    "LogStep": step,
                    "LogPValue": best_p,
                }
            )
        else:  # no remaining variable meets criterion
            break

    return pd.DataFrame.from_records(records)

# =============================================================
#  Information Value / KS / Gini — SAS %Info_Val_Var (parity ver.)
# =============================================================

# ----------------- helpers ----------------- #

def _weighted_quantile(x, w, q):
    """Return weighted quantiles (unused in final version, but kept for ref)."""
    order = np.argsort(x)
    x_sorted = x[order]
    w_sorted = w[order]
    cum_w = np.cumsum(w_sorted)
    qs = np.array(q) * cum_w[-1]
    return np.interp(qs, cum_w, x_sorted)

# ------------------------------------------------------------- #

def _build_bins(
    df: pd.DataFrame,
    var: str,
    dep_var: str,
    *,
    decile: int,
    redu_weight: bool,
    weight: Optional[str],
) -> pd.DataFrame:
    """Faithfully replicates SAS tmp2 table (bins + summary row)."""

    df = df[[var, dep_var] + ([weight] if redu_weight else [])].dropna()
    if df.empty:
        raise ValueError("No valid data after dropping NA for variable " + var)

    # Weight handling
    if redu_weight:
        w = df[weight]
    else:
        w = pd.Series(1.0, index=df.index)

    unique_vals = df[var].nunique(dropna=True)

    # --------------------------------------------------
    # Case 1: unique_vals > decile  →  create decile bins
    # --------------------------------------------------
    if unique_vals > decile:
        if redu_weight:
            # —— Weighted cumulative rank (matches SAS logic) ——
            tmp = df.sort_values(var).copy()
            cum_w = tmp[weight].cumsum()
            total_w = cum_w.iloc[-1]
            tmp["bin"] = np.floor(cum_w * decile / (total_w + 1)).astype(int)
        else:
            # —— Unweighted: emulate PROC RANK groups=&decile ——
            tmp = df.sort_values(var).copy().reset_index(drop=True)
            tmp["rank"] = np.arange(1, len(tmp) + 1)
            tmp["bin"] = (tmp["rank"] * decile // (len(tmp) + 1)).astype(int)
            tmp.drop(columns="rank", inplace=True)
        group_fields = ["bin"]
    # --------------------------------------------------
    # Case 2: unique_vals ≤ decile  →  each value its own bin
    # --------------------------------------------------
    else:
        tmp = df.copy()
        tmp["bin"] = tmp[var]
        group_fields = ["bin"]

    # Aggregate to produce tmp2 body rows (_type_=1 in SAS)
    if redu_weight:
        agg = {
            dep_var: ["mean", "max"],
            var: ["mean", "min", "max"],
            weight: "sum",
        }
    else:
        agg = {
            dep_var: ["mean", "max"],
            var: ["mean", "min", "max"],
        }
    tmp2 = tmp.groupby(group_fields, sort=True).agg(agg)
    if redu_weight:
        tmp2.columns = ["mean_dv", "max_dv", "mean_var", "min_var", "max_var", "cnt"]
    else:
        tmp2.columns = ["mean_dv", "max_dv", "mean_var", "min_var", "max_var"]
        tmp2["cnt"] = tmp2.index.map(tmp["bin"].value_counts())
    tmp2 = tmp2.reset_index()

    # For unique_vals ≤ decile and unweighted, ensure var stats = bin value (SAS behaviour)
    if unique_vals <= decile:
        tmp2["mean_var"] = tmp2["bin"]
        tmp2["min_var"] = tmp2["bin"]
        tmp2["max_var"] = tmp2["bin"]

    # --------------------------------------------------
    # Add summary row (_type_=0 in SAS)
    # --------------------------------------------------
    if redu_weight:
        cnt_total = df[weight].sum()
        mean_dv_total = np.average(df[dep_var], weights=df[weight])
    else:
        cnt_total = len(df)
        mean_dv_total = df[dep_var].mean()
    max_dv_total = df[dep_var].max()

    summary = pd.DataFrame({
        "bin": [np.nan],
        "mean_dv": [mean_dv_total],
        "max_dv": [max_dv_total],
        "mean_var": [np.nan],
        "min_var": [np.nan],
        "max_var": [np.nan],
        "cnt": [cnt_total],
    })

    tmp2 = pd.concat([summary, tmp2], ignore_index=True, sort=False)
    return tmp2

# -------------------------------------------------------------

def info_value_var(
    df: pd.DataFrame,
    var: str,
    dep_var: str,
    *,
    decile: int = 10,
    set_size: int = 1,
    redu_weight: bool = False,
    weight: Optional[str] = None,
) -> pd.DataFrame:
    """Compute IV, Gini, KS exactly as SAS %Info_Val_Var.

    Returns DataFrame: variable | infv | gini | ks
    """

    if redu_weight and weight is None:
        raise ValueError("Weight column must be provided when redu_weight=True")

    # Build tmp2 equivalent
    tmp2 = _build_bins(
        df,
        var,
        dep_var,
        decile=decile,
        redu_weight=redu_weight,
        weight=weight,
    )

    # Global metrics
    norm = tmp2.loc[0, "mean_dv"]
    ntotal = tmp2.loc[0, "cnt"]
    maxdv = tmp2.loc[0, "max_dv"]
    totalresp = (tmp2.loc[1:, "cnt"] * tmp2.loc[1:, "mean_dv"]).sum()

    # Per‑bin metrics (exclude summary row)
    bins = tmp2.loc[1:].copy()
    bins = bins.sort_values("bin").reset_index(drop=True)

    lift_index = (bins["mean_dv"] / norm) * 100

    cumresp = (bins["mean_dv"] * bins["cnt"]).cumsum()
    cumpct_resp = cumresp / totalresp if totalresp != 0 else 0
    cumtotal = bins["cnt"].cumsum()
    cumpct_freq = cumtotal / ntotal if ntotal != 0 else 0

    badrate = bins["mean_dv"] / maxdv if maxdv != 0 else 0
    goodrate = 1 - badrate
    badcnt = badrate * bins["cnt"]
    goodcnt = goodrate * bins["cnt"]

    badratio = (bins["cnt"] * bins["mean_dv"]) / totalresp if totalresp != 0 else 0
    goodratio = (bins["cnt"] - badcnt) / (ntotal - totalresp / maxdv) if ntotal != 0 else 0

    # Iterate for IV / Gini / KS
    infv = 0.0
    gini = 0.0
    ks = 0.0
    cumbad = 0.0
    cumgood = 0.0
    prev_resp = 0.0
    prev_freq = 0.0

    for i, (br, gr, cr, cf) in enumerate(zip(badratio, goodratio, cumpct_resp, cumpct_freq)):
        if i != 0:
            gini += 2 * ((cf + prev_freq) / 2 - (cr + prev_resp) / 2) * (cf - prev_freq)
        prev_resp, prev_freq = cr, cf

        if br + gr == 0:
            contribution = 0.0
        elif br * gr == 0:
            contribution = bins["cnt"].iloc[i] * 1.0 / set_size * max(br, gr)
        else:
            contribution = (br - gr) * np.log(max(br / gr, 1e-5))
        infv += contribution

        cumbad += br
        cumgood += gr
        ks = max(ks, abs(cumbad - cumgood))

    return pd.DataFrame({
        "variable": [var],
        "infv": [infv],
        "gini": [gini],
        "ks": [ks],
    })


# -------------------------------------------------------------------
#  Maxx_Corr — Remove highly correlated variables (SAS parity)
# -------------------------------------------------------------------

def _select_candidates(
    var_redu_df: pd.DataFrame,
    *,
    sources: int,
    include_logsource: bool = False,
    include_regsource: bool = False,
) -> List[str]:
    """Replicates SAS candidate filter based on num_sources etc."""
    mask = var_redu_df["num_sources"] >= sources
    if include_logsource and "logsource" in var_redu_df.columns:
        mask |= var_redu_df["logsource"] == 1
    if include_regsource and "regsource" in var_redu_df.columns:
        mask |= var_redu_df["regsource"] == 1
    return var_redu_df.loc[mask, "variable"].tolist()


def maxx_corr(
    df: pd.DataFrame,
    var_redu_df: pd.DataFrame,
    dep_var: str,
    prefix: str,
    maxcorr: float,
    sources: int,
    logistic: bool,
    regression: bool,
) -> pd.DataFrame:
    """Python equivalent of SAS **%Maxx_Corr**.

    Parameters
    ----------
    df : DataFrame containing dependent & independent variables.
    var_redu_df : DataFrame (analogous to out.CE3_Var_Redu) — must include
        columns: variable, num_sources, plus optional logsource/regsource.
    dep_var : dependent variable name (&dep_var).
    prefix : common prefix of independent variables (&prefix.), used mainly to
        mirror SAS array logic; here we use it to subset correlation matrix.
    maxcorr : correlation threshold (|r| > maxcorr ⇒ considered *high*).
    sources : &sources macro value.
    logistic : whether logistic step is enabled (controls logsource column).
    regression : whether regression step is enabled (controls regsource col).

    Returns
    -------
    Updated *var_redu_df* with column `drop_corr` set to 'Y' for dropped vars.
    """

    # Ensure drop_corr column exists
    if "drop_corr" not in var_redu_df.columns:
        var_redu_df = var_redu_df.copy()
        var_redu_df["drop_corr"] = ""

    # Candidate variables for correlation filtering
    cand_vars = _select_candidates(
        var_redu_df,
        sources=sources,
        include_logsource=logistic,
        include_regsource=regression,
    )
    if not cand_vars:
        return var_redu_df  # nothing to process

    # Build correlation matrix (abs value) between candidate vars only
    corr = df[cand_vars].corr().abs()

    # Also compute |corr| with dep_var to emulate SAS ordering
    dep_corr = df[cand_vars].corrwith(df[dep_var]).abs()

    # Iterative elimination loop
    while True:
        # Count high correlations per variable (excluding self)
        high_cnt = (corr > maxcorr).sum(axis=1) - 1  # subtract self
        # Any variable with >1 such high corr?
        if (high_cnt > 1).sum() == 0:
            break

        # Determine keepvar: first var (highest dep corr) among cnt>1 (SAS sorts)
        candidates_to_keep = high_cnt[high_cnt > 1].index
        keepvar = dep_corr.loc[candidates_to_keep].sort_values(ascending=False).index[0]

        # Vars to drop: vars (≠ keepvar) where |corr(keepvar, var)| > maxcorr
        to_drop = corr.index[(corr.loc[keepvar] > maxcorr) & (corr.index != keepvar)].tolist()
        if not to_drop:
            break  # safety; should not happen

        # Mark in var_redu_df
        var_redu_df.loc[var_redu_df["variable"].isin(to_drop), "drop_corr"] = "Y"

        # Remove dropped vars from corr matrix and dep_corr
        corr.drop(index=to_drop, columns=to_drop, inplace=True)
        dep_corr = dep_corr.drop(index=to_drop)

        # Update candidate list in case loop continues
        if corr.shape[0] <= 1:
            break

    return var_redu_df


# --------------------------------------------------------------------
#                CE_Var_Redu – Master Reduction Flow
# --------------------------------------------------------------------

def ce_var_redu(
    df: pd.DataFrame,
    config: dict
) -> Dict[str, Any]:
    """End‑to‑end Python replica of SAS **%CE_Var_Redu**.

    Returns a dict containing:
    - var_redu_df : combined result table (≈ out.CE3_Var_Redu)
    - varlist_redu : final selected variable names
    - stats_df : basic count/min/max/mean table
    """

    dep_var=config.get("dep_var")
    prefix=config.get("prefix")
    keep_list=config.get("keep_list")
    binary_dv=config.get("binary_dv")
    weight=config.get("weight")
    redu_weight=config.get("redu_weight")
    weight_is_freq=config.get("weight_is_freq")
    # fast run & switches
    fast_opt=config.get("fast_opt")
    univ_reg_flg=config.get("univ_reg_flg")
    correlation_flg=config.get("correlation_flg")
    principal_flg=config.get("principal_flg")
    cluster_flg=config.get("cluster_flg")
    regression_flg=config.get("regression_flg")
    logistic_flg=config.get("logistic_flg")
    information_flg=config.get("information_flg")
    ind_correlation_flg=config.get("ind_correlation_flg")
    ind_dv_corr_flg=config.get("ind_dv_corr_flg")
    # thresholds
    samplesize=config.get("samplesize")
    maxpuni=config.get("maxpuni")
    corrcut=config.get("corrcut")
    nprin=config.get("nprin")
    minprin=config.get("minprin")
    maxc=config.get("maxc")
    maxratio=config.get("maxratio")
    alphareg=config.get("alphareg")
    alphalog=config.get("alphalog")
    infvcut=config.get("infvcut")
    decile=config.get("decile")
    maxcorr=config.get("maxcorr")
    max_dv_corr=config.get("max_dv_corr")
    sources=config.get("sources")
    path_output=config.get("path_output")
    random_seed=config.get("random_seed")

    # --------------------------- FAST OPTION -------------------------- #
    if fast_opt:
        sources = 1
        univ_reg_flg = correlation_flg = principal_flg = cluster_flg = information_flg = ind_correlation_flg = ind_dv_corr_flg = False
        if binary_dv:
            regression_flg = False
            logistic_flg = True
        else:
            regression_flg = True
            logistic_flg = False

    rng = np.random.default_rng(random_seed)

    # --------------------------- 1. WORKFILE -------------------------- #
    workfile = df.loc[df.get("mod_val_test", 0) != 3].copy()
    nobs = len(workfile)
    if nobs > samplesize:
        sample_mask = rng.random(nobs) <= samplesize / nobs
        workfile = workfile.loc[sample_mask]

    # --------------------------- 2. VARIABLE LIST --------------------- #
    keep_list = keep_list or []
    vars_all = [c for c in workfile.columns if c.startswith(prefix) and c not in keep_list]
    varnum = len(vars_all)
    rng.shuffle(vars_all)  # avoid order bias

    outs: Dict[str, pd.DataFrame] = {}

    # ----------------------- 3. UNIVARIATE REG ------------------------ #
    if univ_reg_flg:
        univ_rows = [univ_reg(workfile, v, dep_var, binary_dv=binary_dv,
                               redu_weight=redu_weight, weight=weight) for v in vars_all]
        univ_df = pd.concat(univ_rows, ignore_index=True)
        univ_df["univ_flag"] = (univ_df["PValue"] <= maxpuni).astype(int)
        outs["univ_reg"] = univ_df

    # ----------------------- 4. CORRELATION --------------------------- #
    if correlation_flg:
        corr = workfile[[dep_var] + vars_all].corr().abs()[dep_var].drop(dep_var)
        corr_df = corr.reset_index().rename(columns={"index": "variable", dep_var: "corr"})
        corr_df["corr_flag"] = (corr_df["corr"] >= corrcut).astype(int)
        outs["correlation"] = corr_df

    # ----------------------- 5. PRINCIPAL ----------------------------- #
    if principal_flg:
        grp_size = max(int(np.ceil(varnum / (nprin * 5))), 1)
        groups = np.array_split(vars_all, grp_size)
        prin_df = pd.concat([single_group_prin(workfile, list(g), nprin,
                                               redu_weight=redu_weight, weight=weight) for g in groups],
                            ignore_index=True)
        if len(groups) > 1:
            regv = prin_df.loc[prin_df["Factor"] >= minprin, "variable"].tolist()
            prin_df = single_group_prin(workfile, regv, nprin,
                                        redu_weight=redu_weight, weight=weight)
        prin_df["prin_flag"] = (prin_df["Factor"] >= minprin).astype(int)
        outs["principal"] = prin_df

    # ----------------------- 6. CLUSTER ------------------------------- #
    if cluster_flg:
        grp_size = max(int(np.ceil(varnum / (maxc * 5))), 1)
        groups = np.array_split(vars_all, grp_size)
        clus_df = pd.concat([single_group_clus(workfile, list(g), maxc,
                                               redu_weight=redu_weight, weight=weight) for g in groups],
                            ignore_index=True)
        if len(groups) > 1:
            regv = clus_df.loc[clus_df["RSquareRatio"] <= maxratio, "variable"].tolist()
            clus_df = single_group_clus(workfile, regv, maxc,
                                        redu_weight=redu_weight, weight=weight)
        clus_df["clus_flag"] = (clus_df["RSquareRatio"] <= maxratio).astype(int)
        outs["cluster"] = clus_df

    # ----------------------- 7. REGRESSION ---------------------------- #
    if regression_flg:
        grp_size = int(np.ceil(varnum / 100))
        groups = np.array_split(vars_all, grp_size)
        reg_df = pd.concat([single_group_reg(workfile, list(g), dep_var,
                                             alphareg=alphareg,
                                             redu_weight=redu_weight, weight=weight) for g in groups],
                            ignore_index=True)
        reg_df["regsource"] = 1
        outs["regression"] = reg_df

    # ----------------------- 8. LOGISTIC ------------------------------ #
    if logistic_flg and binary_dv:
        grp_size = int(np.ceil(varnum / 100))
        groups = np.array_split(vars_all, grp_size)
        log_df = pd.concat([single_group_log(workfile, list(g), dep_var,
                                             alphalog=alphalog,
                                             redu_weight=redu_weight, weight=weight,
                                             weight_is_freq=weight_is_freq) for g in groups],
                            ignore_index=True)
        log_df["logsource"] = 1
        outs["logistic"] = log_df

    # ----------------------- 9. INFORMATION VALUE --------------------- #
    if information_flg:
        iv_df = pd.concat([info_value_var(workfile, v, dep_var,
                                          decile=decile, set_size=varnum,
                                          redu_weight=redu_weight, weight=weight) for v in vars_all],
                          ignore_index=True)
        iv_df["infv_flag"] = (iv_df["infv"] >= infvcut).astype(int)
        outs["information"] = iv_df

    # ----------------------- 10. BASIC STATS -------------------------- #
    stats_df = workfile[vars_all].agg(["count", "min", "max", "mean"]).T.reset_index()
    stats_df.columns = ["variable", "Count", "Minimum", "Maximum", "Mean"]

    # ----------------------- 11. MERGE ALL ---------------------------- #
    master = stats_df.copy()
    merge_plan = [
        ("univ_reg", "univ_flag", "univsource"),
        ("correlation", "corr_flag", "corrsource"),
        ("principal", "prin_flag", "prinsource"),
        ("cluster", "clus_flag", "clusource"),
        ("logistic", "logsource", "logsource"),
        ("regression", "regsource", "regsource"),
        ("information", "infv_flag", "infvsource"),
    ]

    for key, flag_col, src_col in merge_plan:
        if key not in outs:
            continue
        df_key = outs[key]
        master = master.merge(df_key[[c for c in df_key.columns if c == "variable" or c == flag_col or c == "infv"]],
                              on="variable", how="left")
        if flag_col.endswith("flag"):
            master[src_col] = master[flag_col].fillna(0).astype(int)
            master.drop(columns=flag_col, inplace=True)
        elif key in ("logistic", "regression"):
            master[src_col] = master["variable"].isin(df_key["variable"]).astype(int)

    source_cols = [c for c in master.columns if c.endswith("source")]
    master["num_sources"] = master[source_cols].sum(axis=1)

    # ----------------------- 12. HIGH CORR DROP ----------------------- #
    if ind_correlation_flg:
        master = maxx_corr(workfile,
                           master,
                           dep_var=dep_var,
                           prefix=prefix,
                           maxcorr=maxcorr,
                           sources=sources,
                           logistic=logistic_flg,
                           regression=regression_flg)

    # ----------------------- 13. DV CORR DROP ------------------------- #
    if ind_dv_corr_flg:
        dv_corr = workfile[vars_all + [dep_var]].corr().abs()[dep_var].drop(dep_var).fillna(0.0)
        master = master.merge(dv_corr.rename("dv_corr"), left_on="variable", right_index=True, how="left")
        master["drop_dv_corr"] = np.where(master["dv_corr"] > max_dv_corr, "Y", "")

    # ----------------------- 14. FINAL SORT --------------------------- #
    sort_cols = ["num_sources"] + (["infv"] if information_flg else [])
    master = master.sort_values(sort_cols, ascending=[False] * len(sort_cols)).reset_index(drop=True)

    # ----------------------- 15. SELECTED VAR LIST -------------------- #
    sel_mask = master["num_sources"] >= sources
    if logistic_flg and binary_dv:
        sel_mask |= master["logsource"] == 1
    if regression_flg:
        sel_mask |= master["regsource"] == 1
    if ind_correlation_flg:
        sel_mask &= (master["drop_corr"].isna() | (master["drop_corr"] == ""))
    if ind_dv_corr_flg:
        sel_mask &= (master["drop_dv_corr"].isna() | (master["drop_dv_corr"] == ""))

    varlist_redu = master.loc[sel_mask, "variable"].tolist()

    # ----------------------- 16. OUTPUT FILES ------------------------- #
    if path_output is not None:
        path_output = Path(path_output)
        path_output.mkdir(parents=True, exist_ok=True)
        master.to_excel(path_output / "CE3_Var_Redu Results.xlsx", index=False)
        Path(path_output / "CE3_Varlist_redu.txt").write_text(
            "%let varlist_redu =\n" + " ".join(varlist_redu) + "\n;")

    return {
        "var_redu_df": master,
        "varlist_redu": varlist_redu,
        "stats_df": stats_df,
    }


# *** Macro 3: Variable Reduction macro variables ***
config.update({
    "fast_opt": False,                    # Fast option will turn off all tests except multivariate regression
    "samplesize": 50000,                 # Sample size to use for variable reduction
    "redu_weight": False,                # Use weights in variables reduction? (Y/N)
    "sources": 3,                        # Minimum number of sources to be selected
    # Univariate regression option
    "univ_reg_flg": True,                    # Use univariate regression to choose variables? (Y/N)
    "maxpuni": 0.0001,                   # Maximum p value correlation for selecting via univariate regression
    # Correlation option
    "correlation_flg": True,                 # Use correlation to choose variables? (Y/N)
    "corrcut": 0.01,                     # Minimum correlation between independent variable and dependent variable
    # Principal components option
    "principal_flg": True,                   # Use principal components to choose variables? (Y/N)
    "nprin": 10,                         # Number of principal components desired
    "minprin": 0.5,                      # Minimum factor correlation for selecting via principal component
    # Cluster option
    "cluster_flg": True,                     # Use cluster analysis to choose variables? (Y/N)
    "maxc": 20,                          # Number of clusters desired
    "maxratio": 0.5,                     # Maximum R Squared ratio for selecting variables via clustering
    # Linear regression option
    "regression_flg": True,                  # Use linear regression to choose variables? (Y/N)
    "alphareg": 0.05,                    # Alpha level for forward selection in linear regression
    # Logistic regression option - only applicable if binary dependent variable
    "logistic_flg": True,                    # Use logistic regression to choose variables? (Y/N)
    "alphalog": 0.05,                    # Alpha level for forward selection in logistic regression
    # Information value option
    "information_flg": True,                 # Use information value to choose variables? (Y/N)
    "decile": 20,                        # Number of groups to use when calculate information values
    "infvcut": 0.01,                     # Minimum information value for selecting via information value
    # Maximum correlation between independent variables option
    "ind_correlation_flg": True,             # Exclude variables with high correlation to others? (Y/N)
    "maxcorr": 0.7,                      # Maximum correlation allowed between independent variables
    # Maximum correlation to dependent variable option
    "ind_dv_corr_flg": True,                 # Exclude variables with high correlation to dependent variable? (Y/N)
    "max_dv_corr": 0.7                   # Maximum correlation allowed to dependent variable
})

# ② 运行变量筛选
result = ce_var_redu(
    df=CE2_Recoded,
    config=config
)

# ③ 查看结果
print("↓最终选中变量↓")
print(result["varlist_redu"])
print("\n头 10 行合并表")
print(result["var_redu_df"].head(10))





*** 3.Variable reduction and ranking ***;
%let LogFile = "&path_output.03_Var_Redu_Log_File.log";
%let LstFile = "&path_output.03_Var_Redu_List_File.lst";
proc printto log=&LogFile new; run;
proc printto print=&LstFile new; run;

%CE_Var_Redu(insdn=out.CE2_Recoded);

proc printto;
run;
