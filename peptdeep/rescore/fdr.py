import numba
import numpy as np
import pandas as pd

@numba.njit
def fdr_to_q_values(
    fdr_values:np.ndarray
)->np.ndarray:
    """convert FDR values to q_values.

    Parameters
    ----------
    fdr_values : np.ndarray
        FDR values, they should be 
        sorted according to the descending order of the `score`

    Returns
    -------
    np.ndarray
        q_values

    """
    q_values = np.zeros_like(fdr_values)
    min_q_value = np.max(fdr_values)
    for i in range(len(fdr_values) - 1, -1, -1):
        fdr = fdr_values[i]
        if fdr < min_q_value:
            min_q_value = fdr
        q_values[i] = min_q_value
    return q_values

def calc_fdr(
    df:pd.DataFrame, 
    score_column:str, 
    decoy_column:str='decoy'
)->pd.DataFrame:
    """Calculate FDR values (q_values in fact) for the given dataframe

    Parameters
    ----------
    df : pd.DataFrame
        PSM dataframe to calculate FDRs

    score_column : str
        score column to sort in decending order

    decoy_column : str, optional
        decoy column in the dataframe. 
        1=target, 0=decoy. Defaults to 'decoy'.

    Returns
    -------
    pd.DataFrame
        PSM dataframe with 'fdr' column added

    """
    df = df.reset_index(drop=True).sort_values(
        [score_column,decoy_column], ascending=False
    )
    target_values = 1-df[decoy_column].values
    decoy_cumsum = np.cumsum(df[decoy_column].values)
    target_cumsum = np.cumsum(target_values)
    fdr_values = decoy_cumsum/target_cumsum
    df['fdr'] = fdr_to_q_values(fdr_values)
    return df

#wrapper
calc_fdr_for_df = calc_fdr

@numba.njit
def fdr_from_ref(
    sorted_scores:np.ndarray, 
    ref_scores:np.ndarray, 
    ref_fdr_values:np.ndarray
)->np.ndarray:
    """ Calculate FDR values from the given reference scores and fdr_values. 
    It is used to extend peptide-level or sequence-level FDR (reference) 
    to each PSM, as PSMs are more useful for quantification.

    Parameters
    ----------
    sorted_scores : np.array
        the scores to calculate FDRs, 
        they must be sorted in decending order.

    ref_scores : np.array
        reference scores that used to 
        calculate ref_fdr_values, also sorted in decending order.

    ref_fdr_values : np.array
        fdr values corresponding to ref_scores

    Returns
    -------
    np.array
        fdr values corresponding to sorted_scores.

    """
    q_values = np.zeros_like(sorted_scores)
    i,j = 0,0
    while i < len(sorted_scores) and j < len(ref_scores):
        if sorted_scores[i] >= ref_scores[j]:
            q_values[i] = ref_fdr_values[j]
            i += 1
        else:
            j += 1
    while i < len(sorted_scores):
        q_values[i] = ref_fdr_values[-1]
        i += 1
    return q_values

def calc_fdr_from_ref(
    df: pd.DataFrame,
    ref_scores:np.ndarray, 
    ref_fdr_values:np.ndarray,
    score_column:str, 
    decoy_column:str='decoy'
)->pd.DataFrame:
    """ Calculate FDR values for a PSM dataframe from the given reference
     scores and fdr_values. It is used to extend peptide-level or 
     sequence-level FDR (reference) to each PSM, as PSMs are more useful 
     for quantification.
    ``

    Parameters
    ----------
    df : pd.DataFrame
        PSM dataframe

    ref_scores : np.array
        reference scores that used to 
        calculate ref_fdr_values, also sorted in decending order.

    ref_fdr_values : np.array
        fdr values corresponding to ref_scores

    score_column : str
        score column in the dataframe

    decoy_column : str, optional
        decoy column in the dataframe. 
        1=target, 0=decoy. Defaults to 'decoy'.

    Returns
    -------
    pd.DataFrame
        dataframe with 'fdr' column added

    """
    df = df.reset_index(drop=True).sort_values(
        [score_column,decoy_column], ascending=False
    )
    sorted_idxes = np.argsort(ref_fdr_values)
    ref_scores = ref_scores[sorted_idxes]
    ref_q_values = ref_fdr_values[sorted_idxes]
    df['fdr'] = fdr_from_ref(
        df.score.values, ref_scores, ref_q_values
    )
    return df

calc_fdr_from_ref_for_df = calc_fdr_from_ref
