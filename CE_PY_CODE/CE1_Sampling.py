import pandas as pd
import numpy as np
import logging
import os
from typing import Tuple, Optional

__all__ = ["CE_Sampling"]


def sample_option_1(
    inds: pd.DataFrame, samplesize: int, seed: Optional[int] = None
) -> pd.DataFrame:

    return inds.sample(n=samplesize, replace=False, random_state=seed)


def sample_option_2(
    inds: pd.DataFrame, samplesize: int, seed: Optional[int] = None
) -> pd.DataFrame:

    orig = inds.copy()
    extra = samplesize - len(orig)
    if extra > 0:
        extra_df = inds.sample(n=extra, replace=True, random_state=seed)
        return pd.concat([orig, extra_df], ignore_index=True)
    return orig


def sample_core(
    inds: pd.DataFrame,
    response_var: str,
    samplesize: int,
    ratio: float,
    min_num_resp: int,
    seed: Optional[int] = None,
    logger: Optional[logging.Logger] = None,
) -> pd.DataFrame:

    # 1. 计算最小 non-response
    min_nonresp = int((min_num_resp / ratio) - min_num_resp)
    if logger:
        logger.debug(f"min_num_nonresp = {min_nonresp}")

    # 2. 计算 n1, n2 并边界校正
    n1 = int(samplesize * ratio)
    n2 = samplesize - n1
    if n1 < min_num_resp:
        if logger:
            logger.warning(f"Too few responses; set n1 = {min_num_resp}")
        n1 = min_num_resp
    if n2 < min_nonresp:
        if logger:
            logger.warning(f"Too few non-responses; set n2 = {min_nonresp}")
        n2 = min_nonresp

    # 3. 分组 responders / non-responders
    df_resp = inds[inds[response_var] > 0]
    df_non = inds[inds[response_var] <= 0]
    if logger:
        logger.debug(
            f"Original responders: {len(df_resp)}, non-responders: {len(df_non)}"
        )

    # 4. 采样或复制
    if len(df_resp) == n1:
        _out1 = df_resp.copy()
    elif len(df_resp) > n1:
        _out1 = sample_option_1(df_resp, n1, seed)
    else:
        _out1 = sample_option_2(df_resp, n1, seed)

    if len(df_non) == n2:
        _out2 = df_non.copy()
    elif len(df_non) > n2:
        _out2 = sample_option_1(df_non, n2, seed)
    else:
        _out2 = sample_option_2(df_non, n2, seed)

    # 5. 合并输出
    result = pd.concat([_out1, _out2], ignore_index=True)
    if logger:
        logger.debug(f"Sampled total rows = {len(result)}")
    return result


def CE_Sampling(
    config: dict, logger: Optional[logging.Logger] = None
) -> Tuple[pd.DataFrame, pd.DataFrame]:

    outds_name = str(config.get("outds_name", "CE1_Resampled"))
    if logger:
        logger.info(f"=== Start CE_Sampling (outds={outds_name}) ===")
    # General Macro Variables
    path_output = config["path_output"]
    inds = config["inds"]
    id = config["id"]
    dep_var = config["dep_var"]
    binary_dv = config["binary_dv"].upper() == "Y"
    weight = config["weight"]
    # Macro 1: Sampling macro variables
    split_portion = config["split_portion"]
    exclusion_if = config["exclusion_if"]
    DS_present = (
        config["DS_present"].upper() == "Y"
    )
    # Optional Macro Variables: These all have defaults that can be used
    # General Macro Variables
    prefix = config.get("prefix", "R1_")
    keep_list = config.get("keep_list", [])
    # Macro 1: Sampling macro variables
    path_DS = config.get("path_DS", "/mnt/projects/shared/pst_qmgisi/Modeling/CE/")
    bootstrap = config.get("bootstrap", "N").upper() == "Y"
    oversampled_rr = config.get("oversampled_rr", 0.1)
    min_num_resp = config.get("min_num_resp", 2000)
    seed = config.get("seed", 123456)

    # 1. Exclusion
    mask_excl = eval(exclusion_if)
    if logger:
        logger.info(f"Excluding {mask_excl.sum()} rows")
    tmp = inds.loc[~mask_excl]

    # 3. Split into _mod, _test and set mod_val_test
    np.random.seed(seed)
    mask_split = np.random.rand(len(tmp)) < split_portion
    mod = tmp.loc[mask_split].reset_index(drop=True)
    mod = mod.assign(mod_val_test=1)
    test = tmp.loc[~mask_split].reset_index(drop=True)
    test = test.assign(mod_val_test=3)
    if logger:
        logger.info(f"Split: model={len(mod)}, test={len(test)}")

    # 4. Bootstrap sampling on _mod
    if binary_dv and bootstrap:
        total = len(mod)
        true_resp = mod[dep_var].mean()
        num_resp = total * true_resp
        if num_resp >= min_num_resp:
            spsize = int((total * true_resp) / oversampled_rr)
        else:
            spsize = int(min_num_resp / oversampled_rr)
            if logger:
                logger.warning(f"Low resp; use min_num_resp={min_num_resp}")
        if logger:
            logger.info(f"Bootstrap spsize={spsize}, rr={oversampled_rr}")
        mod = sample_core(
            mod, dep_var, spsize, oversampled_rr, min_num_resp, seed, logger
        )
        if logger:
            logger.info(
                f"Resampled model freq: {mod[dep_var].value_counts().to_dict()}"
            )

    # 5. Combine and Validation carve-out
    outds = pd.concat([mod, test], ignore_index=True)
    rng2 = np.random.RandomState(seed + 1 if seed is not None else None)
    idx_mod = outds[outds["mod_val_test"] == 1].index
    # SAS 按 50% ranuni 随机分配 Val
    val_mask = rng2.rand(len(idx_mod)) > 0.5
    val_idx = idx_mod[val_mask]
    outds.loc[val_idx, "mod_val_test"] = 2
    if logger:
        logger.info(f"Marked {len(val_idx)} rows as Validation")

    # 6. Compute sample-rate table
    rate = (
        outds.groupby("mod_val_test")[dep_var]
        .mean()
        .reset_index(name="rate")
        .assign(
            partition=lambda df_: df_["mod_val_test"].map(
                {1: "Model", 2: "Val", 3: "Test"}
            )
        )
        .loc[:, ["partition", "rate"]]
    )

    # 7. Write outputs
    os.makedirs(path_output, exist_ok=True)
    outds.to_csv(os.path.join(path_output, f"{outds_name}.csv"), index=False)
    rate.to_csv(os.path.join(path_output, f"{outds_name}_Sample_Rate.csv"), index=False)
    if logger:
        logger.info(
            f"Wrote {outds.shape[0]} rows to {outds_name}.csv and sample-rate table"
        )
        logger.info(f"=== End CE_Sampling (outds={outds_name}) ===")

    return outds, rate
