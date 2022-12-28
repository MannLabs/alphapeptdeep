import numpy as np
import numba
import pandas as pd
import tqdm

from alphabase.peptide.fragment import (
    create_fragment_mz_dataframe, 
    get_charged_frag_types
)
from peptdeep.mass_spec.ms_reader import (
    ms2_reader_provider, MSReaderBase
)

@numba.njit
def match_centroid_mz(
    spec_mzs:np.ndarray, 
    query_mzs:np.ndarray, 
    spec_mz_tols:np.ndarray,
)->np.ndarray:
    """
    Matching query masses against sorted MS2/spec centroid masses, 
    only closest (minimal abs mass error) peaks are returned.
    
    Parameters
    ----------
    spec_mzs : np.ndarray
        MS2 or spec mz values, 1-D float array

    query_mzs : np.ndarray
        query mz values, n-D float array

    spec_mz_tols : np.ndarray
        Da tolerance array, same shape as spec_mzs

    Returns
    -------
    np.ndarray
        np.ndarray of int32, the shape is the same as query_mzs.
        -1 means no peaks are matched for the query mz

    """
    idxes = np.searchsorted(spec_mzs, query_mzs)
    ret_indices = np.empty_like(query_mzs, dtype=np.int32)
    # ret_indices[:] = -1
    for i,idx in np.ndenumerate(idxes):
        min_merr = abs(spec_mzs[idx-1]-query_mzs[i])
        min_idx = -1
        if min_merr <= spec_mz_tols[idx-1]:
            min_idx = idx-1
        if idx < len(spec_mzs):
            merr = abs(spec_mzs[idx]-query_mzs[i])
            if merr <= spec_mz_tols[idx] and merr < min_merr: 
                min_idx = idx
        ret_indices[i] = min_idx
    return ret_indices


@numba.njit
def match_profile_mz(
    spec_mzs:np.ndarray,
    query_mzs:np.ndarray,
    spec_mz_tols:np.ndarray,
    spec_intens:np.ndarray,
)->np.ndarray:
    """
    Matching query masses against sorted MS2/spec profile masses,
    only highest peaks are returned.

    Parameters
    ----------
    spec_mzs : np.ndarray
        MS2 or spec mz values, 1-D float array

    query_mzs : np.ndarray
        query mz values, n-D float array

    spec_mz_tols : np.ndarray
        Da tolerance array, same shape as spec_mzs

    Returns
    -------
    np.ndarray
        np.ndarray of int32, the shape is the same as query_mzs.
        -1 means no peaks are matched for the query mz

    """
    idxes = np.searchsorted(spec_mzs, query_mzs)
    ret_indices = np.empty_like(query_mzs, dtype=np.int32)
    for i,idx in np.ndenumerate(idxes):
        if idx == len(spec_mzs):
            idx = idx-1

        highest = 0
        highest_idx = -1
            
        for _idx in range(idx, -1, -1):
            if abs(spec_mzs[_idx]-query_mzs[i])>spec_mz_tols[_idx]:
                break
            if highest < spec_intens[_idx]:
                highest = spec_intens[_idx]
                highest_idx = _idx
        for _idx in range(idx+1, len(spec_mzs)):
            if abs(spec_mzs[_idx]-query_mzs[i])>spec_mz_tols[_idx]:
                break
            if highest < spec_intens[_idx]:
                highest = spec_intens[_idx]
                highest_idx = _idx
        ret_indices[i] = highest_idx
    return ret_indices

@numba.njit
def match_first_last_profile_mz(
    spec_mzs:np.ndarray,
    query_mzs:np.ndarray,
    spec_mz_tols:np.ndarray,
)->np.ndarray:
    """
    Matching query masses against sorted MS2/spec profile masses,
    both first and last m/z values are returned.

    Parameters
    ----------
    spec_mzs : np.ndarray
        MS2 or spec mz values, 1-D float array

    query_mzs : np.ndarray
        query mz values, n-D float array

    spec_mz_tols : np.ndarray
        Da tolerance array, same shape as spec_mzs

    Returns
    -------
    np.ndarray
        2D np.ndarray of int32 with first and last matched index 
        for the query mz. The shape is the same as (len(query_mzs),2).
        -1 means no peaks are matched for the query mz
    """
    idxes = np.searchsorted(spec_mzs, query_mzs)
    first_indices = np.empty_like(
        query_mzs, dtype=np.int32
    )
    last_indices = np.empty_like(
        query_mzs, dtype=np.int32
    )
    first_indices[:] = -1
    last_indices[:] = -1
    for i,idx in np.ndenumerate(idxes):
        if idx == len(spec_mzs):
            idx = idx-1
        for _idx in range(idx, -1, -1):
            if spec_mzs[_idx]<query_mzs[i]-spec_mz_tols[_idx]:
                break
            else:
                first_indices[i] = _idx
        for _idx in range(idx, len(spec_mzs)):
            if abs(spec_mzs[_idx]-query_mzs[i])>spec_mz_tols[_idx]:
                break
            else:
                last_indices[i] = _idx
    return first_indices, last_indices


@numba.njit
def match_one_raw_with_numba(
    spec_idxes, frag_start_idxes, frag_stop_idxes,
    all_frag_mzs,
    all_spec_mzs, all_spec_intensities, 
    peak_start_idxes, peak_end_idxes,
    matched_intensities, matched_mz_errs,
    ppm, tol,
):
    """ 
    Internel function to match fragment mz values to spectrum mz values.
    Matched_mz_errs[i] = np.inf if no peaks are matched.
    """
    for spec_idx, frag_start, frag_end in zip(
        spec_idxes, frag_start_idxes, frag_stop_idxes
    ):
        peak_start = peak_start_idxes[spec_idx]
        peak_end = peak_end_idxes[spec_idx]
        if peak_end == peak_start: continue
        spec_mzs = all_spec_mzs[peak_start:peak_end]
        spec_intens = all_spec_intensities[peak_start:peak_end]
        
        if ppm:
            spec_mz_tols = spec_mzs*tol*1e-6
        else:
            spec_mz_tols = np.full_like(spec_mzs, tol)

        frag_mzs = all_frag_mzs[frag_start:frag_end,:].copy()
        
        matched_idxes = match_centroid_mz(
            spec_mzs, frag_mzs, spec_mz_tols
        ).reshape(-1)

        matched_intens = spec_intens[matched_idxes]
        matched_intens[matched_idxes==-1] = 0

        matched_mass_errs = np.abs(
            spec_mzs[
                matched_idxes.reshape(-1)
            ]-frag_mzs.reshape(-1)
        )
        matched_mass_errs[matched_idxes==-1] = np.inf

        matched_intensities[
            frag_start:frag_end,:
        ] = matched_intens.reshape(frag_mzs.shape)

        matched_mz_errs[
            frag_start:frag_end,:
        ] = matched_mass_errs.reshape(frag_mzs.shape)


class PepSpecMatch(object):
    """Main entry for peptide-spectrum matching"""
    def __init__(self,
        charged_frag_types = get_charged_frag_types(
            ['b','y','b_modloss','y_modloss'],
            2
        )
    ):
        self.charged_frag_types = charged_frag_types

    def _preprocess_psms(self, psm_df):
        pass

    def get_fragment_mz_df(self, psm_df):
        return create_fragment_mz_dataframe(
            psm_df, self.charged_frag_types
        )

    def match_ms2_one_raw(self, 
        psm_df_one_raw: pd.DataFrame,
        ms2_file:str,
        ms2_file_type:str='alphapept',
        ppm:bool=True, tol:float=20.0,
    )->tuple:
        """Matching psm_df_one_raw against ms2_file

        Parameters
        ----------
        psm_df_one_raw : pd.DataFrame
            psm dataframe 
            that contains only one raw file

        ms2_file : str
            ms2 file path

        ms2_file_type : str, optional
            ms2 file type, could be ["thermo","alphapept","mgf"].
            Default to 'alphapept'

        ppm : bool, optional
            if use ppm tolerance. Defaults to True.

        tol : float, optional
            tolerance value. Defaults to 20.0.

        Returns
        -------
        tuple:
            pd.DataFrame: psm dataframe with fragment index information.
            
            pd.DataFrame: fragment mz dataframe.
            
            pd.DataFrame: matched intensity dataframe.
            
            pd.DataFrame: matched mass error dataframe. 
            np.inf if a fragment is not matched.
            
        """
        self._preprocess_psms(psm_df_one_raw)
        psm_df = psm_df_one_raw
        if isinstance(ms2_file, MSReaderBase):
            ms2_reader = ms2_file
        else:
            ms2_reader = ms2_reader_provider.get_reader(
                ms2_file_type
            )
            ms2_reader.load(ms2_file)

        add_spec_info_list = []
        if 'rt_norm' not in psm_df.columns:
            add_spec_info_list.append('rt')

        if (
            'mobility' not in psm_df.columns and 
            'mobility' in ms2_reader.spectrum_df.columns
        ):
            add_spec_info_list.append('mobility')

        if len(add_spec_info_list) > 0:
            # pfind does not report RT in the result file
            psm_df = psm_df.reset_index().merge(
                ms2_reader.spectrum_df[['spec_idx']+add_spec_info_list],
                how='left',
                on='spec_idx',
            ).set_index('index')

            if 'rt' in add_spec_info_list:
                psm_df['rt_norm'] = psm_df.rt/ms2_reader.spectrum_df.rt.max()

        fragment_mz_df = self.get_fragment_mz_df(psm_df)
        
        matched_intensity_df = pd.DataFrame(
            np.zeros_like(
                fragment_mz_df.values, dtype=np.float64
            ), 
            columns=fragment_mz_df.columns
        )

        matched_mz_err_df = pd.DataFrame(
            np.full_like(
                fragment_mz_df.values, np.inf, dtype=np.float64
            ), 
            columns=fragment_mz_df.columns
        )
        
        for (
            spec_idx, frag_start_idx, frag_stop_idx
        ) in psm_df[[
            'spec_idx', 'frag_start_idx', 
            'frag_stop_idx'
        ]].values:
            (
                spec_mzs, spec_intens
            ) = ms2_reader.get_peaks(spec_idx)
            if len(spec_mzs)==0: continue

            if ppm:
                mz_tols = spec_mzs*tol*1e-6
            else:
                mz_tols = np.full_like(spec_mzs, tol)

            frag_mzs = fragment_mz_df.values[
                frag_start_idx:frag_stop_idx,:
            ]
            
            matched_idxes = match_centroid_mz(
                spec_mzs, frag_mzs, mz_tols
            )
            matched_intens = spec_intens[matched_idxes]
            matched_intens[matched_idxes==-1] = 0

            matched_mz_errs = np.abs(
                spec_mzs[matched_idxes]-frag_mzs
            )
            matched_mz_errs[matched_idxes==-1] = np.inf

            matched_intensity_df.values[
                frag_start_idx:frag_stop_idx,:
            ] = matched_intens

            matched_mz_err_df.values[
                frag_start_idx:frag_stop_idx,:
            ] = matched_mz_errs

        return (
            psm_df, fragment_mz_df, 
            matched_intensity_df, matched_mz_err_df
        )

    def _match_ms2_centroid_one_raw(self, raw_name, df_group):
        if raw_name in self._ms2_file_dict:
            if isinstance(self._ms2_file_dict[raw_name], MSReaderBase):
                ms2_reader = self._ms2_file_dict[raw_name]
            else:
                ms2_reader = ms2_reader_provider.get_reader(
                    self._ms2_file_type
                )
                ms2_reader.load(self._ms2_file_dict[raw_name])
            if self.rt_not_in_df:
                # pfind does not report RT in the result file
                _df = df_group.reset_index().merge(
                    ms2_reader.spectrum_df[['spec_idx','rt']],
                    how='left',
                    on='spec_idx',
                ).set_index('index')

                _df['rt_norm'] = _df.rt/ms2_reader.spectrum_df.rt.max()
                self.psm_df.loc[
                    _df.index, ['rt','rt_norm']
                ] = _df[['rt','rt_norm']]

            match_one_raw_with_numba(
                df_group.spec_idx.values,
                df_group.frag_start_idx.values,
                df_group.frag_stop_idx.values,
                self.fragment_mz_df.values,
                ms2_reader.peak_df.mz.values, 
                ms2_reader.peak_df.intensity.values,
                ms2_reader.spectrum_df.peak_start_idx.values,
                ms2_reader.spectrum_df.peak_end_idx.values,
                self.matched_intensity_df.values,
                self.matched_mz_err_df.values,
                self.ppm, self.tol
            )
    
    def match_ms2_centroid(self,
        psm_df: pd.DataFrame,
        ms2_file_dict: dict, #raw_name: ms2_file_path or ms_reader object
        ms2_file_type:str = 'alphapept', # or 'mgf', or 'thermo'
        ppm=True, tol=20.0,
    ):
        """Matching PSM dataframe against the ms2 files in ms2_file_dict
        This method will store matched values as attributes:
        - self.psm_df
        - self.fragment_mz_df
        - self.matched_intensity_df
        - self.matched_mz_err_df

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM dataframe

        ms2_file_dict : dict
            {raw_name: ms2 path}

        ms2_file_type : str, optional
            Could be 'alphapept', 'mgf' or 'thermo'. 
            Defaults to 'alphapept'.

        ppm : bool, optional
            Defaults to True.

        tol : float, optional
            PPM units, defaults to 20.0.
            
        """
        self._preprocess_psms(psm_df)
        self.psm_df = psm_df

        if 'frag_start_idx' in self.psm_df.columns:
            del self.psm_df['frag_start_idx']
            del self.psm_df['frag_stop_idx']
            
        self.fragment_mz_df = self.get_fragment_mz_df(self.psm_df)
        
        self.matched_intensity_df = pd.DataFrame(
            np.zeros_like(
                self.fragment_mz_df.values, dtype=np.float64
            ), 
            columns=self.fragment_mz_df.columns
        )

        self.matched_mz_err_df = pd.DataFrame(
            np.full_like(
                self.fragment_mz_df.values, np.inf, dtype=np.float64
            ), 
            columns=self.fragment_mz_df.columns
        )
        
        self._ms2_file_dict = ms2_file_dict
        self._ms2_file_type = ms2_file_type
        self.ppm = ppm
        self.tol = tol

        if 'rt_norm' not in self.psm_df.columns:
            self.rt_not_in_df = True
        else:
            self.rt_not_in_df = False
        for raw_name, df_group in tqdm.tqdm(
            self.psm_df.groupby('raw_name')
        ):
            self._match_ms2_centroid_one_raw(raw_name, df_group)
