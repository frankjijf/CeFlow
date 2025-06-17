import pandas as pd
import numpy as np
import logging
import os

__all__ = ['CE_Sampling']

def _setup_logger(log_path: str, lst_path: str = None):
    """Configure a logger to write to log and list files at specified paths."""
    logger = logging.getLogger('CE_Sampling')
    logger.setLevel(logging.DEBUG)
    # Clear existing handlers
    for h in logger.handlers[:]:
        logger.removeHandler(h)
    # Ensure directory exists for log
    os.makedirs(os.path.dirname(log_path), exist_ok=True)
    fh = logging.FileHandler(log_path, mode='w')
    fh.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s %(levelname)s: %(message)s')
    fh.setFormatter(formatter)
    logger.addHandler(fh)
    # Optional list file handler
    if lst_path:
        os.makedirs(os.path.dirname(lst_path), exist_ok=True)
        lh = logging.FileHandler(lst_path, mode='w')
        lh.setLevel(logging.INFO)
        lh.setFormatter(formatter)
        logger.addHandler(lh)
    return logger


def sample_option_1(df: pd.DataFrame, n: int, seed: int = None) -> pd.DataFrame:
    """Replicates SAS Sample_Option_1: simple random sample without replacement."""
    return df.sample(n=n, replace=False, random_state=seed)


def sample_option_2(df: pd.DataFrame, n: int, seed: int = None) -> pd.DataFrame:
    """Replicates SAS Sample_Option_2: if desired > available, bootstrap extra."""
    orig = df.copy()
    extra = n - len(df)
    if extra > 0:
        extra_df = df.sample(n=extra, replace=True, random_state=seed)
        return pd.concat([orig, extra_df], ignore_index=True)
    return orig


def sample_core(df: pd.DataFrame,
                response_var: str,
                samplesize: int,
                ratio: float,
                min_num_resp: int,
                seed: int = None,
                logger: logging.Logger = None) -> pd.DataFrame:
    """
    Full replication of SAS %Sample_Core including warnings and stratified sampling.
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
            logger.warning(f"Too few non-responses; setting non-responses to {min_nonresp}")
        n2 = min_nonresp

    df_resp = df[df[response_var] > 0]
    df_non  = df[df[response_var] <= 0]
    if logger:
        logger.debug(f"Original responders: {len(df_resp)}, non-responders: {len(df_non)}")

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


def CE_Sampling(df: pd.DataFrame,
                dep_var: str,
                split_frac: float = 0.7,
                exclusion_fn = None,
                DS_present: bool = False,
                ds_recodes_fn = None,
                binary_dv: bool = False,
                bootstrap: bool = False,
                oversampled_rr: float = 0.05,
                min_num_resp: int = 500,
                seed: int = None,
                path_output: str = './',
                log_file: str = 'CE_Sampling.log',
                lst_file: str = 'CE_Sampling.lst') -> (pd.DataFrame, pd.DataFrame):
    """
    Apple-to-apple replication of SAS %CE_Sampling with proper log and list file paths.
    Returns resampled DataFrame and sample rate summary DataFrame.
    """
    # Ensure output directory exists
    os.makedirs(path_output, exist_ok=True)
    # Build full log/list paths
    log_path = os.path.join(path_output, log_file)
    lst_path = os.path.join(path_output, lst_file) if lst_file else None
    # Setup logging with correct paths
    logger = _setup_logger(log_path, lst_path)
    logger.info("Starting CE_Sampling")

    df_work = df.copy()
    # Exclusion
    if exclusion_fn is not None:
        mask_excl = exclusion_fn(df_work)
        logger.info(f"Excluding {mask_excl.sum()} records")
        df_work = df_work.loc[~mask_excl]
    # DS Recodes
    if DS_present and ds_recodes_fn:
        df_work = ds_recodes_fn(df_work)
        logger.info("Applied DS recodes")

    # Split into mod/test
    rng_split = np.random.RandomState(seed)
    mask_mod = rng_split.rand(len(df_work)) < split_frac
    df_work['mod_val_test'] = np.where(mask_mod, 1, 3)
    mod = df_work[df_work['mod_val_test'] == 1].reset_index(drop=True)
    test = df_work[df_work['mod_val_test'] == 3].reset_index(drop=True)
    logger.info(f"Split: {len(mod)} model, {len(test)} test")

    # Original freq QC
    freq_orig = df_work[dep_var].value_counts()
    logger.info(f"Original {dep_var} frequencies: {freq_orig.to_dict()}")

    # Bootstrap oversampling
    if binary_dv and bootstrap:
        temp_count = mod[dep_var].count()
        true_resp = mod[dep_var].mean()
        num_resp = temp_count * true_resp
        if num_resp >= min_num_resp:
            spsize = int((temp_count * true_resp) / oversampled_rr)
        else:
            spsize = int(min_num_resp / oversampled_rr)
            logger.warning(f"Too few responses; using min_num_resp={min_num_resp}")
        logger.info(f"spsize={spsize}, oversampled_rr={oversampled_rr}, true_resp={true_resp}")
        mod = sample_core(mod, dep_var, spsize, oversampled_rr, min_num_resp, seed, logger)
        freq_mod = mod[dep_var].value_counts()
        freq_test = test[dep_var].value_counts()
        logger.info(f"Training oversampled {dep_var} frequencies: {freq_mod.to_dict()}")
        logger.info(f"Test {dep_var} frequencies: {freq_test.to_dict()}")

    # Combine and mark validation
    res = pd.concat([mod, test], ignore_index=True)
    rng_val = np.random.RandomState(seed+1 if seed is not None else None)
    mask_val = (res['mod_val_test'] == 1) & (rng_val.rand(len(res)) > 0.5)
    res.loc[mask_val, 'mod_val_test'] = 2
    logger.info(f"Marked {mask_val.sum()} records as Validation")

    # Sample rate summary
    rate = res.groupby('mod_val_test')[dep_var].mean().reset_index(name='rate')
    mapping = {1: 'Model', 2: 'Val', 3: 'Test'}
    rate['partition'] = rate['mod_val_test'].map(mapping)
    sample_rate_df = rate[['partition', 'rate']]
    # Write summary CSV
    sr_path = os.path.join(path_output, 'CE1_Sample_Rate.csv')
    sample_rate_df.to_csv(sr_path, index=False)
    logger.info(f"Wrote sample rate summary to {sr_path}")
    logger.info("Completed CE_Sampling")
    return res, sample_rate_df

# Load Titanic dataset
df_titanic = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/Titanic-Dataset.csv')

# Run CE_Sampling on Titanic (no bootstrap)
resampled_df, sample_rate_df = CE_Sampling(
    df=df_titanic,
    dep_var='Survived',
    split_frac=0.7,
    binary_dv=True,
    bootstrap=False,
    seed=654321,
    path_output='D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/'
)

resampled_df.to_csv(
    'D:/OneDrive/CE_PROJECT/Python_CE_Project_GH/SAS2PYTHON/CE1_Resampled_Titanic.csv',
    index=False
)