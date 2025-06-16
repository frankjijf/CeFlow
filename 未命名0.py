# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 14:09:47 2025

@author: frank
"""

import pandas as pd
import numpy as np
from typing import Optional, Union


def sample_without_replacement(df: pd.DataFrame, n: int, seed: Optional[int] = None) -> pd.DataFrame:
    """Sample n rows without replacement."""
    return df.sample(n=n, replace=False, random_state=seed)


def sample_with_replacement(df: pd.DataFrame, n: int, seed: Optional[int] = None) -> pd.DataFrame:
    """Sample n rows with replacement, ensuring original df included once."""
    original = df.copy()
    extra_n = n - len(original)
    if extra_n <= 0:
        return original
    rng = np.random.RandomState(seed)
    extra = original.sample(n=extra_n, replace=True, random_state=rng)
    return pd.concat([original, extra], ignore_index=True)


def sample_core(
    df: pd.DataFrame,
    response_var: str,
    samplesize: int,
    ratio: float,
    seed: Optional[int] = None,
    min_num_resp: Optional[int] = None
) -> pd.DataFrame:
    """
    Core sampling replicating SAS Sample_Core macro:
      - ratio: desired response rate in sample (oversampled_rr)
      - samplesize: total desired sample size (spsize)
      - min_num_resp: minimum number of responders
    """
    # Desired counts
    n_resp = int(samplesize * ratio)
    n_nonresp = samplesize - n_resp

    # Enforce minimum responder and non-responder counts if provided
    if min_num_resp is not None:
        min_nonresp = int(min_num_resp / ratio - min_num_resp)
        if n_resp < min_num_resp:
            n_resp = min_num_resp
            print(f"Warning: Too few responders, setting to min_num_resp={min_num_resp}")
        if n_nonresp < min_nonresp:
            n_nonresp = min_nonresp
            print(f"Warning: Too few non-responders, setting to min_nonresp={min_nonresp}")

    # Partition
    df_resp = df[df[response_var] > 0]
    df_non = df[df[response_var] <= 0]

    # Sample responders
    if len(df_resp) >= n_resp:
        out_resp = sample_without_replacement(df_resp, n_resp, seed)
    else:
        out_resp = sample_with_replacement(df_resp, n_resp, seed)

    # Sample non-responders
    if len(df_non) >= n_nonresp:
        out_non = sample_without_replacement(df_non, n_nonresp, seed)
    else:
        out_non = sample_with_replacement(df_non, n_nonresp, seed)

    return pd.concat([out_resp, out_non], ignore_index=True)


def ce_sampling(
    df: pd.DataFrame,
    dep_var: str,
    binary_dv: bool = True,
    split_frac: float = 0.7,
    exclusion_condition: Optional[str] = None,
    ds_present: bool = False,
    path_ds: Optional[str] = None,
    bootstrap: bool = False,
    oversampled_rr: float = 0.05,
    min_num_resp: int = 500,
    seed_split: int = 123456,
    seed_final: int = 9545
) -> pd.DataFrame:
    """
    Replicates the CE_Sampling macro from SAS:
      1. Exclude rows based on exclusion_condition
      2. Initial split into modeling (1) and test (3) by split_frac
      3. Optionally apply DS recodes if ds_present=True
      4. If binary_dv & bootstrap: oversample modeling set to achieve oversampled_rr
      5. Combine modeling & test, then split modeling into model (1) and validation (2)
      6. Returns DataFrame with column 'mod_val_test' (1=model, 2=validation, 3=test)
    """
    df = df.copy()

    # 1. Exclusion
    if exclusion_condition:
        df = df.query(f"not ({exclusion_condition})")

    # 2. Initial split
    rng = np.random.RandomState(seed_split)
    df['mod_val_test'] = np.where(rng.rand(len(df)) < split_frac, 1, 3)

    # 3. DS recodes stub
    if ds_present:
        # TODO: implement DS recodes loading from path_ds
        pass

    # Separate modeling and test
    mod = df[df['mod_val_test'] == 1].copy()
    test = df[df['mod_val_test'] == 3].copy()

    # 4. Bootstrap/oversampling if requested
    if binary_dv and bootstrap:
        # Compute original response rate
        true_resp = mod[dep_var].mean()
        count1 = mod.shape[0]
        num_resp = int(count1 * true_resp)
        if num_resp >= min_num_resp:
            spsize = int(num_resp / oversampled_rr)
        else:
            spsize = int(min_num_resp / oversampled_rr)
            print(f"Warning: Too few responses. Number of Responses set to {min_num_resp}")
        # Perform core sampling
        mod = sample_core(
            mod,
            response_var=dep_var,
            samplesize=spsize,
            ratio=oversampled_rr,
            seed=seed_split,
            min_num_resp=min_num_resp
        )

    # 5. Combine and split into model/validation/test
    out = pd.concat([mod, test], ignore_index=True)
    rng2 = np.random.RandomState(seed_final)
    mask = (out['mod_val_test'] == 1) & (rng2.rand(len(out)) > 0.5)
    out.loc[mask, 'mod_val_test'] = 2

    return out


# Load and apply
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_250616/Titanic-Dataset.csv')
sampled = ce_sampling(df, dep_var='Survived', binary_dv=True, split_frac=0.7, bootstrap=False)

from ace_tools import display_dataframe_to_user
display_dataframe_to_user('Sampled Titanic Data (first 10 rows)', sampled.head(10))

partition_counts = sampled['mod_val_test'].value_counts().rename_axis('partition').reset_index(name='count')
display_dataframe_to_user('Partition Counts', partition_counts)