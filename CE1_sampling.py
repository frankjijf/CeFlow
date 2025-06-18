"""
CE1_sampling.py

This script replicates SAS sampling methods in Python, including simple random sampling, stratified sampling, and oversampling techniques. It is designed to handle binary dependent variables and supports logging and exclusion logic.

Functions:
- _setup_logger: Configures logging to specified log and list files.
- sample_option_1: Implements simple random sampling without replacement.
- sample_option_2: Handles oversampling when the desired sample size exceeds available records.
- sample_core: Performs stratified sampling with warnings for insufficient responses/non-responses.
- CE_Sampling: Main function replicating SAS %CE_Sampling, including exclusion logic, DS recodes, and partitioning into model, validation, and test sets.

Usage:
- Load a dataset (e.g., Titanic-Dataset.csv).
- Call CE_Sampling with appropriate parameters.
- Outputs include resampled data and sample rate summary.

"""

import pandas as pd
import numpy as np
import logging
import os
from typing import Tuple, Optional
import re

__all__ = ["CE_Sampling"]


def _setup_logger(log_path: str, lst_path: Optional[str] = None):
    """
    Configure a logger to write to log and list files at specified paths.

    Parameters:
    - log_path (str): Path to the log file.
    - lst_path (str, optional): Path to the list file. Defaults to None.

    Returns:
    - logging.Logger: Configured logger instance.
    """
    logger = logging.getLogger("CE_Sampling")
    logger.setLevel(logging.DEBUG)
    # Clear existing handlers
    for h in logger.handlers[:]:
        logger.removeHandler(h)
    # Ensure directory exists for log
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    fh = logging.FileHandler(log_path, mode="w")
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter("%(asctime)s %(levelname)s: %(message)s")
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    # Optional list file handler
    if lst_path:
        os.makedirs(os.path.dirname(lst_path), exist_ok=True)
        lh = logging.FileHandler(lst_path, mode="w")
        lh.setLevel(logging.INFO)
        lh.setFormatter(formatter)
        logger.addHandler(lh)
    return logger


def sample_option_1(df: pd.DataFrame, n: int, seed: Optional[int] = None) -> pd.DataFrame:
    """
    Implements simple random sampling without replacement.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - n (int): Number of samples to draw.
    - seed (int, optional): Random seed for reproducibility. Defaults to None.

    Returns:
    - pd.DataFrame: Sampled DataFrame.
    """
    return df.sample(n=n, replace=False, random_state=seed)


def sample_option_2(df: pd.DataFrame, n: int, seed: Optional[int] = None) -> pd.DataFrame:
    """
    Handles oversampling when the desired sample size exceeds available records.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - n (int): Desired sample size.
    - seed (int, optional): Random seed for reproducibility. Defaults to None.

    Returns:
    - pd.DataFrame: Oversampled DataFrame.
    """
    orig = df.copy()
    extra = n - len(df)
    if extra > 0:
        extra_df = df.sample(n=extra, replace=True, random_state=seed)
        return pd.concat([orig, extra_df], ignore_index=True)
    return orig


def sample_core(
    df: pd.DataFrame,
    response_var: str,
    samplesize: int,
    ratio: float,
    min_num_resp: int,
    seed: Optional[int] = None,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:
    """
    Performs stratified sampling with warnings for insufficient responses/non-responses.

    Parameters:
    - df (pd.DataFrame): Input DataFrame.
    - response_var (str): Column name of the response variable.
    - samplesize (int): Total sample size.
    - ratio (float): Desired ratio of responses to non-responses.
    - min_num_resp (int): Minimum number of responses required.
    - seed (int, optional): Random seed for reproducibility. Defaults to None.
    - logger (logging.Logger, optional): Logger instance for logging. Defaults to None.

    Returns:
    - pd.DataFrame: Stratified sampled DataFrame.
    """
    min_nonresp = int((min_num_resp / ratio) - min_num_resp)
    if logger:
        logger.debug(f"min_num_nonresp set to {min_nonresp}")

    n1 = int(samplesize * ratio)
    n2 = samplesize - n1
    if n1 < min_num_resp:
        if logger:
            logger.warning(f"Too few responses; setting responses to {min_num_resp}")
        n1 = min_num_resp
    if n2 < min_nonresp:
        if logger:
            logger.warning(
                f"Too few non-responses; setting non-responses to {min_nonresp}"
            )
        n2 = min_nonresp

    df_resp = df[df[response_var] > 0]
    df_non = df[df[response_var] <= 0]
    if logger:
        logger.debug(
            f"Original responders: {len(df_resp)}, non-responders: {len(df_non)}"
        )

    if len(df_resp) == n1:
        out1 = df_resp.copy()
    elif len(df_resp) > n1:
        out1 = sample_option_1(df_resp, n1, seed)
    else:
        out1 = sample_option_2(df_resp, n1, seed)

    if len(df_non) == n2:
        out2 = df_non.copy()
    elif len(df_non) > n2:
        out2 = sample_option_1(df_non, n2, seed)
    else:
        out2 = sample_option_2(df_non, n2, seed)

    result = pd.concat([out1, out2], ignore_index=True)
    if logger:
        logger.debug(f"Sampled total records: {len(result)}")
    return result


def CE_Sampling(df: pd.DataFrame, config: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    # 1. 解析 config
    dep_var      = config["dep_var"]
    seed         = config.get("seed", None)
    binary_dv    = str(config.get("binary_dv","N")).upper() == "Y"
    bootstrap    = str(config.get("bootstrap","N")).upper() == "Y"
    oversampled_rr = float(config.get("oversampled_rr", 0.05))
    min_num_resp = int(config.get("min_num_resp", 500))
    # split_if 解析阈值
    split_if_str = config.get("split_if", "np.random.rand() < 0.7")
    m = re.search(r"<\s*([0-9]*\.?[0-9]+)", split_if_str)
    split_frac = float(m.group(1)) if m else 0.7
    exclusion_if_str = config.get("exclusion_if", "")
    DS_present   = str(config.get("DS_Present","N")).upper() == "Y"
    ds_recodes_fn= config.get("ds_recodes_fn", None)
    path_output  = config.get("path_output", "./")
    log_file     = config.get("log_file", "CE_Sampling.log")
    lst_file     = config.get("lst_file", "CE_Sampling.lst")
    val_frac     = float(config.get("val_frac", 0.5))
    
    # 2. Logger
    log_path = os.path.join(path_output, log_file)
    lst_path = os.path.join(path_output, lst_file) if lst_file else None
    logger   = _setup_logger(log_path, lst_path)
    logger.info("=== Start CE_Sampling ===")

    df_work = df.copy()

    # 3. 排除逻辑
    if exclusion_if_str:
        mask_excl = df_work.eval(exclusion_if_str)
        logger.info(f"Excluding {mask_excl.sum()} rows by exclusion_if")
        df_work = df_work.loc[~mask_excl]

    # 4. DS recodes
    if DS_present and ds_recodes_fn:
        df_work = ds_recodes_fn(df_work)
        logger.info("Applied DS recodes")

    # 5. 划分 mod/test
    rng = np.random.RandomState(seed)
    mask_mod = rng.rand(len(df_work)) < split_frac
    df_work["mod_val_test"] = np.where(mask_mod, 1, 3)
    mod  = df_work[df_work["mod_val_test"]==1].reset_index(drop=True)
    test = df_work[df_work["mod_val_test"]==3].reset_index(drop=True)
    logger.info(f"Split: {len(mod)} model, {len(test)} test")

    # 原始频次
    freq_orig = df_work[dep_var].value_counts()
    logger.info(f"Original {dep_var} freq: {freq_orig.to_dict()}")


    # 6. Bootstrap oversampling
    if binary_dv and bootstrap:
        total = len(mod)
        true_resp = mod[dep_var].mean()
        num_resp  = total * true_resp
        if num_resp >= min_num_resp:
            spsize = int((total * true_resp) / oversampled_rr)
        else:
            spsize = int(min_num_resp / oversampled_rr)
            logger.warning(f"Low resp; use min_num_resp={min_num_resp}")
        logger.info(f"Bootstrap spsize={spsize}, rr={oversampled_rr}, true_resp={true_resp}")
        mod = sample_core(mod, dep_var, spsize, oversampled_rr, min_num_resp, seed, logger)
        logger.info(f"Oversampled model freq: {mod[dep_var].value_counts().to_dict()}")
        logger.info(f"Test freq: {test[dep_var].value_counts().to_dict()}")

    # 7. 合并并打标 Validation
    res = pd.concat([mod, test], ignore_index=True)
    rng2 = np.random.RandomState(seed+1 if seed is not None else None)
    idx_mod = res.loc[res["mod_val_test"]==1].index
    n_val   = int(len(idx_mod)*val_frac)
    val_idx = rng2.choice(idx_mod, size=n_val, replace=False) if n_val>0 else []
    res.loc[val_idx, "mod_val_test"] = 2
    logger.info(f"Marked {len(val_idx)} as Validation (val_frac={val_frac})")

    # 8. 样本率汇总
    rate = (
        res.groupby("mod_val_test")[dep_var]
           .mean()
           .reset_index(name="rate")
           .assign(partition = lambda df_: df_["mod_val_test"].map({1:"Model",2:"Val",3:"Test"}))
           .loc[:, ["partition","rate"]]
    )
    out_csv = os.path.join(path_output, "CE1_Sample_Rate.csv")
    rate.to_csv(out_csv, index=False)
    logger.info(f"Wrote sample-rate to {out_csv}")

    logger.info("=== End CE_Sampling ===")
    return res, rate