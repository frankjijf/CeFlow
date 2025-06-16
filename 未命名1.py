# -*- coding: utf-8 -*-
"""
Created on Mon Jun 16 14:50:25 2025

@author: frank
"""

import pandas as pd
import numpy as np
from typing import Dict, Tuple
from scipy import stats


def pnum(
    df: pd.DataFrame,
    var: str,
#    dep_var: str,
    impmethod: str = 'median',
    cap_floor: bool = True,
    transform: bool = True
) -> Dict[str, float]:
    """
    Mimics SAS pnum macro: computes var_LB, var_UB, var_median, var_miss (imputation value).
    """
    series = df[var].dropna()
    p1, p25, p75, p99 = np.percentile(series, [1, 25, 75, 99])
    var_min, var_max = series.min(), series.max()
    iqr = max(p75 - p25, p99 - p75, p25 - p1)
    var_lb = min(max(p25 - 1.5 * iqr, var_min), p1)
    var_ub = max(min(p75 + 1.5 * iqr, var_max), p99)
    if var_lb == var_ub:
        var_lb, var_ub = var_min, var_max
    var_mid = (var_max + var_min) / 2
    var_median = series.median()
    var_mean = series.mean()
    # Determine var_miss based on impmethod
    m = impmethod.upper()
    if m in ('MEAN', 'STD'):
        var_miss = var_mean
    elif m in ('MEDIAN', 'IQR', 'MAD'):
        var_miss = var_median
    elif m == 'RANGE':
        var_miss = var_min
    elif m == 'MIDRANGE':
        var_miss = var_mid
    elif m in ('SUM', 'EUCLEN', 'USTD', 'MAXABS'):
        var_miss = 0.0
    else:
        # fallback to median
        var_miss = var_median
    return {
        'var_LB': var_lb,
        'var_UB': var_ub,
        'var_median': var_median,
        'var_miss': var_miss
    }


def prof1(
    df: pd.DataFrame,
    var: str,
    dep_var: str
) -> pd.DataFrame:
    """
    Mimics SAS prof1: for ordinal variables with few categories.
    Returns DataFrame with columns: variable, xcategory, xcount, xmean
    """
    grp = (
        df.groupby(var)[dep_var]
          .agg(xcount='size', xmean='mean')
          .reset_index()
    )
    grp['xcategory'] = grp[var].astype(str).fillna('Missing')
    grp['variable'] = var
    return grp[['variable', 'xcategory', 'xcount', 'xmean']]


def prof2(
    df: pd.DataFrame,
    var: str,
    dep_var: str,
    num_category: int = 10,
    equal_dist: bool = False
) -> pd.DataFrame:
    """
    Mimics SAS prof2: bins continuous/ordinal into num_category groups
    Returns DataFrame with lo, hi, xcount, xmean for each bin
    """
    temp = df[[var, dep_var]].copy()
    if equal_dist:
        bins = np.linspace(temp[var].min(), temp[var].max(), num_category+1)
        temp['bin'] = pd.cut(temp[var], bins=bins, include_lowest=True)
    else:
        temp['bin'] = pd.qcut(temp[var], q=num_category, duplicates='drop')
    prof = (
        temp.dropna(subset=['bin'])
            .groupby('bin')
            .agg(lo=(var, 'min'), hi=(var, 'max'), xcount=(var, 'size'), xmean=(dep_var, 'mean'))
            .reset_index()
    )
    prof['variable'] = var
    # label xcategory
    def label_row(i,row):
        if row['bin'].isna(): return 'Missing'
        if isinstance(row['bin'], pd.Interval):
            return f"{row['bin'].left:.2f} to {row['bin'].right:.2f}"
        return str(row['bin'])
    prof['xcategory'] = prof.apply(lambda row: label_row(None,row), axis=1)
    return prof[['variable', 'xcategory', 'xcount', 'xmean', 'lo', 'hi']]


def prof3(
    prof_df: pd.DataFrame,
    var: str,
    overall_avg: float,
    nobs: int,
    typical_star: Tuple[float, float] = (110, 90)
) -> pd.DataFrame:
    """
    Mimics SAS prof3: combine profiling results into final profile table
    Adds percent, index, star.
    """
    df = prof_df.copy()
    df['percent'] = df['xcount'] / nobs
    df['Average_DV'] = df['xmean']
    df['index'] = df['Average_DV'] / overall_avg * 100
    high, low = typical_star
    def star(idx):
        if idx >= high: return '* (+)'
        if idx > 100: return '  (+)'
        if idx <= low: return '* (-)'
        return '  (-)'
    df['star'] = df['index'].apply(star)
    df['variable'] = var
    return df[['variable', 'xcategory', 'xcount', 'percent', 'Average_DV', 'index', 'star']]


# Load test data
df = pd.read_csv('D:/OneDrive/CE_PROJECT/Python_CE_Project_250616/Titanic-Dataset.csv')

# Test pnum on 'Age'
pnum_stats = pnum(df, 'Age')
print("pnum stats for Age:", pnum_stats)

# Test prof1 on 'Pclass'
prof1_df = prof1(df, 'Pclass', 'Survived')
from ace_tools import display_dataframe_to_user
display_dataframe_to_user('prof1 (Pclass)', prof1_df)

# Test prof2 on 'Age'
prof2_df = prof2(df, 'Age', 'Survived', num_category=5, equal_dist=False)
display_dataframe_to_user('prof2 (Age, 5 bins)', prof2_df)

# Test prof3 on result of prof2
overall_avg = df['Survived'].mean()
nobs = len(df)
prof3_df = prof3(prof2_df, 'Age', overall_avg, nobs)
display_dataframe_to_user('prof3 (Age)', prof3_df)