"""Utilities for psm matching. Moved from alpharaw."""

from typing import Union

import numba
import numpy as np

from alpharaw.match.match_utils import match_closest_peaks, match_highest_peaks
from alpharaw.ms_data_base import MSData_Base, ms_reader_provider


@numba.jit(nogil=True)
def match_one_raw_with_numba(
    spec_idxes: np.ndarray,
    frag_start_idxes: np.ndarray,
    frag_stop_idxes: np.ndarray,
    all_frag_mzs: np.ndarray,
    all_frag_mz_tols: np.ndarray,
    all_spec_mzs: np.ndarray,
    all_spec_intensities: np.ndarray,
    peak_start_idxes: np.ndarray,
    peak_stop_idxes: np.ndarray,
    matched_intensities: np.ndarray,
    matched_mz_errs: np.ndarray,
    match_closest: bool = True,
) -> None:
    """
    Internel function to match fragment mz values to spectrum mz values.
    Matched_mz_errs[i] = np.inf if no peaks are matched.

    Results will saved in place of matched_intensities
    and matched_mz_errs.
    """
    for spec_idx, frag_start, frag_stop in zip(
        spec_idxes, frag_start_idxes, frag_stop_idxes
    ):
        if spec_idx == -1:
            continue
        peak_start = peak_start_idxes[spec_idx]
        peak_stop = peak_stop_idxes[spec_idx]
        if peak_stop == peak_start:
            continue
        spec_mzs = all_spec_mzs[peak_start:peak_stop]
        spec_intens = all_spec_intensities[peak_start:peak_stop]

        frag_mzs = all_frag_mzs[frag_start:frag_stop, :].copy()
        frag_mz_tols = all_frag_mz_tols[frag_start:frag_stop, :].copy()

        if match_closest:
            matched_idxes = match_closest_peaks(
                spec_mzs, spec_intens, frag_mzs, frag_mz_tols
            ).reshape(-1)
        else:
            matched_idxes = match_highest_peaks(
                spec_mzs, spec_intens, frag_mzs, frag_mz_tols
            ).reshape(-1)

        matched_intens = spec_intens[matched_idxes]
        matched_intens[matched_idxes == -1] = 0

        matched_mass_errs = np.abs(
            spec_mzs[matched_idxes.reshape(-1)] - frag_mzs.reshape(-1)
        )
        matched_mass_errs[matched_idxes == -1] = np.inf

        matched_intensities[frag_start:frag_stop, :] = matched_intens.reshape(
            frag_mzs.shape
        )

        matched_mz_errs[frag_start:frag_stop, :] = matched_mass_errs.reshape(
            frag_mzs.shape
        )


def load_ms_data(
    ms_file: Union[str, MSData_Base],
    ms_file_type: str = "alpharaw_hdf",
    process_count: int = 8,
) -> MSData_Base:
    """
    Load MS file and get `MSData_Base` object.

    Parameters
    ----------
    ms_file : str | MSData_Base
        ms2 file path

    ms_file_type : str, optional
        ms2 file type, can be
        ["alpharaw_hdf", "thermo", "sciex", "alphapept_hdf", "mgf"].
        Default to 'alpharaw_hdf'.

    Returns
    -------
    MSData_Base:
        Instance of sub-class of `MSData_Base`.
    """
    if isinstance(ms_file, MSData_Base):
        return ms_file
    else:
        raw_data = ms_reader_provider.get_reader(
            ms_file_type, process_count=process_count
        )
        if raw_data is None:
            print(
                f"Raw file type '{ms_file_type}' is not registered in 'ms_reader_provider'"
            )
            return
        raw_data.import_raw(ms_file)
        return raw_data


@numba.njit
def get_ion_count_scores(
    frag_mz_values: np.ndarray,
    matched_intens: np.ndarray,
    frag_start_idxes: np.ndarray,
    frag_stop_idxes: np.ndarray,
    min_mz: float = 200,
):
    """
    TODO Deprecated
    """
    scores = []
    for i in range(len(frag_start_idxes)):
        scores.append(
            np.count_nonzero(
                matched_intens[frag_start_idxes[i] : frag_stop_idxes[i], :]
                .copy()
                .reshape(-1)[
                    frag_mz_values[frag_start_idxes[i] : frag_stop_idxes[i], :]
                    .copy()
                    .reshape(-1)
                    >= min_mz
                ]
            )
        )
    return np.array(scores, np.int32)
