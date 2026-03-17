"""This module provides the `PepSpecMatch` class to match fragment ions of PSMs. Moved from alpharaw."""

from typing import Tuple, Union

import numba
import numpy as np
import pandas as pd
import tqdm
from alphabase.peptide.fragment import (
    create_fragment_mz_dataframe,
    get_charged_frag_types,
)

from peptdeep.match.utils import load_ms_data, match_one_raw_with_numba
from alpharaw import register_all_readers
from normal_dia import NormalDIAGrouper
from alpharaw.match.match_utils import (
    match_closest_peaks,
    match_highest_peaks,
)
from alpharaw.match.spec_finder import find_dia_spec_idxes_same_window
from alpharaw.ms_data_base import (
    PEAK_INTENSITY_DTYPE,
    PEAK_MZ_DTYPE,
    MSData_Base,
)
from alpharaw.utils.ms_path_utils import parse_ms_files_to_dict

register_all_readers()  # TODO remove this import side effect


class PepSpecMatch:
    """
    Extract fragment ions from MS2 data.
    The extracted information can be used for visualization of peak annotation or
    PeptDeep transfer learnining for the MS2 model.
    """

    match_closest: bool = True
    use_ppm: bool = True
    #: matching mass tolerance
    tolerance: float = 20.0
    ms_loader_thread_num: int = 4

    def __init__(
        self,
        charged_frag_types: list = None,
        match_closest: bool = True,
        use_ppm: bool = True,
        tol_value: float = 20.0,
    ):
        """
        Parameters
        ----------
        charged_frag_types : list, optional
            fragment types with charge states,
            e.g. ['b_z1', 'y_z2', 'b_modloss_z1', 'y_H2O_z2'].
            Defaults to `get_charged_frag_types(['b','y','b_modloss','y_modloss'], 2)`.

        match_closest : bool, optional
            if True, match the closest peak for a m/z;
            if False, matched the higest peak for a m/z in the tolerance range.
            By default True.

        use_ppm : bool, optional
            If use ppm other wise Da, by default True.

        tol_value : float, optional
            Matching tolerance value (ppm or Da based on `use_ppm`)
            for peak annotation, by default 20.0
        """
        self.charged_frag_types = (
            get_charged_frag_types(["b", "y", "b_modloss", "y_modloss"], 2)
            if charged_frag_types is None
            else charged_frag_types
        )
        self.match_closest = match_closest
        self.use_ppm = use_ppm
        self.tolerance = tol_value

    def get_fragment_mz_df(self) -> pd.DataFrame:
        """
        Call :func:`alphabase.peptide.fragment.create_fragment_mz_dataframe`
        for :attr:`PepSpecMatch.psm_df` and :attr:`PepSpecMatch.charged_frag_types`.


        Returns
        -------
        DataFrame
            The fragment m/z dataframe in alphabase format.
        """
        return create_fragment_mz_dataframe(
            self.psm_df,
            self.charged_frag_types,
            dtype=PEAK_MZ_DTYPE,
        )

    def _add_missing_columns_to_psm_df(
        self, psm_df: pd.DataFrame, raw_data: MSData_Base = None
    ):
        """
        Add missing "rt", "nce", "rt_norm", ("mobility") columns to `psm_df` inplace if missing.

        Parameters
        ----------
        psm_df : pd.DataFrame
            psm dataframe to be processed.
        raw_data : MSData_Base, optional
            The `MSData_Base`. If None, `self.raw_data`. by default None.

        Returns
        -------
        DataFrame
            The original `psm_df` with missing columns added.
        """
        if raw_data is None:
            raw_data = self.raw_data
        add_spec_info_list = []
        if "rt" not in psm_df.columns:
            add_spec_info_list.append("rt")

        if "nce" in raw_data.spectrum_df.columns and "nce" not in psm_df.columns:
            add_spec_info_list.append("nce")

        if (
            "mobility" not in psm_df.columns
            and "mobility" in raw_data.spectrum_df.columns
        ):
            add_spec_info_list.append("mobility")

        if len(add_spec_info_list) > 0:
            # pfind does not report RT in the result file
            psm_df = (
                psm_df.reset_index()
                .merge(
                    raw_data.spectrum_df[["spec_idx"] + add_spec_info_list],
                    how="left",
                    on="spec_idx",
                )
                .set_index("index")
            )

            if "rt" in add_spec_info_list:
                psm_df["rt_norm"] = psm_df.rt / raw_data.spectrum_df.rt.max()
        # if 'rt_sec' not in psm_df.columns:
        #     psm_df['rt_sec'] = psm_df.rt*60
        return psm_df

    def _prepare_matching_dfs(self) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Prepare empty `fragment_mz_df`, `matched_intensity_df`,
        and `matched_mz_err_df` dataframes to extract peak matching information
        for `self.psm_df`. These three dataframes will be only used internally
        in :class:`PepSpecMatch` objects.

        Returns
        -------
        Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame]
            pd.DataFrame: fragment mz dataframe.

            pd.DataFrame: intensity dataframe to match.

            pd.DataFrame: mz error dataframe to match.
        """
        fragment_mz_df = self.get_fragment_mz_df()

        matched_intensity_df = pd.DataFrame(
            np.zeros_like(fragment_mz_df.values, dtype=PEAK_INTENSITY_DTYPE),
            columns=fragment_mz_df.columns,
        )

        matched_mz_err_df = pd.DataFrame(
            np.full_like(fragment_mz_df.values, np.inf, dtype=PEAK_MZ_DTYPE),
            columns=fragment_mz_df.columns,
        )
        return (fragment_mz_df, matched_intensity_df, matched_mz_err_df)

    def load_ms_data(
        self,
        ms_file: Union[str, MSData_Base],
        ms_file_type: str = "alpharaw_hdf",
        process_count: int = 8,
        **kwargs,
    ):
        """Load MS file to set `self.raw_data`.

        Parameters
        ----------
        ms_file : str | MSData_Base
            Absolute or relative path of the ms2 file.

        ms_file_type : str, optional
            ms2 file type, could be
            ["alpharaw_hdf","thermo","sciex","alphapept_hdf","mgf"].
            Default to 'alpharaw_hdf'.
        """
        self.raw_data = load_ms_data(ms_file, ms_file_type, process_count=process_count)

    def get_peaks(self, spec_idx: int, **kwargs):
        return self.raw_data.get_peaks(spec_idx)

    def _match_one_psm(
        self,
        peak_mzs: np.ndarray,
        peak_intensities: np.ndarray,
        fragment_mz_df: pd.DataFrame,
        matched_intensity_df: pd.DataFrame,
        matched_mz_err_df: pd.DataFrame,
        frag_start_idx: int,
        frag_stop_idx: int,
    ):
        """
        Match fragments of one precursor (located by `frag_start_idx` and `frag_stop_idx`)
        against the corresponding `peak_mzs`.

        Parameters
        ----------
        peak_mzs : np.ndarray
            Peak m/z values to be matched.
        peak_intensities : np.ndarray
            Peak intensities to be matched.
        fragment_mz_df : pd.DataFrame
            fragment m/z dataframe to be matched.
        matched_intensity_df : pd.DataFrame
            The dataframe to store matched intensity values.
        matched_mz_err_df : pd.DataFrame
            The dataframe to store matched mz error values.
        frag_start_idx : int
            fragment start index of the given PSM.
        frag_stop_idx : int
            fragment stop index of the given PSM.
        """
        if len(peak_mzs) == 0:
            return

        peak_mzs = peak_mzs.astype(PEAK_MZ_DTYPE)

        frag_mzs = fragment_mz_df.values[frag_start_idx:frag_stop_idx, :]

        if self.use_ppm:
            mz_tols = frag_mzs * self.tolerance * 1e-6
        else:
            mz_tols = np.full_like(frag_mzs, self.tolerance)

        if self.match_closest:
            matched_idxes = match_closest_peaks(
                peak_mzs, peak_intensities, frag_mzs, mz_tols
            )
        else:
            matched_idxes = match_highest_peaks(
                peak_mzs,
                peak_intensities,
                frag_mzs,
                mz_tols,
            )

        matched_intens = peak_intensities[matched_idxes]
        matched_intens[matched_idxes == -1] = 0

        matched_mz_errs = np.abs(peak_mzs[matched_idxes] - frag_mzs)
        matched_mz_errs[matched_idxes == -1] = np.inf

        matched_intensity_df.values[frag_start_idx:frag_stop_idx, :] = matched_intens

        matched_mz_err_df.values[frag_start_idx:frag_stop_idx, :] = matched_mz_errs

    def match_ms2_one_raw(
        self,
        psm_df_one_raw: pd.DataFrame,
        verbose: bool = False,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Matching psm_df_one_raw against self.raw_data
        after `self.load_ms_data()`

        Parameters
        ----------
        psm_df_one_raw : pd.DataFrame
            psm dataframe
            that contains only one raw file

        Returns
        -------
        Tuple:
            pd.DataFrame: psm dataframe with fragment index information.

            pd.DataFrame: fragment mz dataframe.

            pd.DataFrame: matched intensity dataframe.

            pd.DataFrame: matched mass error dataframe.
              np.inf if a fragment is not matched.
        """
        self.psm_df = psm_df_one_raw

        psm_df_one_raw = self._add_missing_columns_to_psm_df(
            psm_df_one_raw, self.raw_data
        )

        (
            fragment_mz_df,
            matched_intensity_df,
            matched_mz_err_df,
        ) = self._prepare_matching_dfs()

        psm_iters = psm_df_one_raw[
            ["spec_idx", "frag_start_idx", "frag_stop_idx"]
        ].values
        if verbose:
            psm_iters = tqdm.tqdm(psm_iters)

        for spec_idx, frag_start_idx, frag_stop_idx in psm_iters:
            (spec_mzs, spec_intens) = self.get_peaks(spec_idx)

            self._match_one_psm(
                spec_mzs,
                spec_intens,
                fragment_mz_df,
                matched_intensity_df,
                matched_mz_err_df,
                frag_start_idx,
                frag_stop_idx,
            )

        return (psm_df_one_raw, fragment_mz_df, matched_intensity_df, matched_mz_err_df)

    def _match_ms2_one_raw_numba(self, raw_name: str, psm_df_one_raw: pd.DataFrame):
        if raw_name in self._ms_file_dict:
            raw_data = load_ms_data(
                self._ms_file_dict[raw_name],
                self._ms_file_type,
                process_count=self.ms_loader_thread_num,
            )
            if raw_data is None:
                return

            psm_df_one_raw = self._add_missing_columns_to_psm_df(
                psm_df_one_raw, raw_data
            )

            if self.use_ppm:
                all_frag_mz_tols = self.fragment_mz_df.values * self.tolerance * 1e-6
            else:
                all_frag_mz_tols = np.full_like(
                    self.fragment_mz_df.values, self.tolerance
                )

            match_one_raw_with_numba(
                psm_df_one_raw.spec_idx.values,
                psm_df_one_raw.frag_start_idx.values,
                psm_df_one_raw.frag_stop_idx.values,
                self.fragment_mz_df.values,
                all_frag_mz_tols,
                raw_data.peak_df.mz.values,
                raw_data.peak_df.intensity.values,
                raw_data.spectrum_df.peak_start_idx.values,
                raw_data.spectrum_df.peak_stop_idx.values,
                self.matched_intensity_df.values,
                self.matched_mz_err_df.values,
                self.match_closest,
            )
        else:
            print(f"`{raw_name}` is not found in ms_file_dict.")
            return
        return psm_df_one_raw

    def match_ms2_multi_raw(
        self,
        psm_df: pd.DataFrame,
        ms_files: Union[dict, list],
        ms_file_type: str = "alpharaw_hdf",
        process_num: int = 1,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Matching PSM dataframe against the ms2 files in ms_files
        This method will store matched values as attributes:
        - self.psm_df
        - self.fragment_mz_df
        - self.matched_intensity_df
        - self.matched_mz_err_df

        Parameters
        ----------
        psm_df : pd.DataFrame
            PSM dataframe

        ms_files : dict | list
            if dict: {raw_name: ms2 path}
            if list: [ms2 path1, ms2 path2]

        ms_file_type : str, optional
            Could be 'alpharaw_hdf', 'mgf' or 'thermo', 'sciex', 'alphapept_hdf'.
            Defaults to 'alphapept'.

        Returns
        -------
        Tuple:
            pd.DataFrame: psm dataframe with fragment index information.

            pd.DataFrame: fragment mz dataframe.

            pd.DataFrame: matched intensity dataframe.

            pd.DataFrame: matched mass error dataframe.
                np.inf if a fragment is not matched.

        """
        self.psm_df = psm_df

        (
            self.fragment_mz_df,
            self.matched_intensity_df,
            self.matched_mz_err_df,
        ) = self._prepare_matching_dfs()

        if isinstance(ms_files, dict):
            self._ms_file_dict = ms_files
        else:
            self._ms_file_dict = parse_ms_files_to_dict(ms_files)

        self._ms_file_type = ms_file_type

        self.ms_loader_thread_num = process_num
        # if process_num <= 1 or len(self._ms_file_dict) <= 1:
        psm_df_list = []
        for raw_name, df_group in tqdm.tqdm(self.psm_df.groupby("raw_name")):
            _df = self._match_ms2_one_raw_numba(raw_name, df_group)
            if _df is not None:
                psm_df_list.append(_df)
        # else:
        #     with mp.get_context("spawn").Pool(processes=process_num) as p:
        #         df_groupby = self.psm_df.groupby('raw_name')
        #         def gen_group_df(df_groupby):
        #             for raw_name, df_group in df_groupby:
        #                 yield (raw_name, df_group)
        #         process_bar(
        #             p.imap_unordered(
        #                 self._match_ms2_one_raw_numba,
        #                 gen_group_df(df_groupby)
        #             ), df_groupby.ngroups
        #         )

        self.psm_df = pd.concat(psm_df_list, ignore_index=True)
        return (
            self.psm_df,
            self.fragment_mz_df,
            self.matched_intensity_df,
            self.matched_mz_err_df,
        )


class PepSpecMatch_DIA(PepSpecMatch):
    """
    Peak annotation for DIA data.
    """

    max_spec_per_query: int = 3
    min_frag_mz: float = 200.0

    def _add_missing_columns_to_psm_df(self, psm_df: pd.DataFrame, raw_data=None):
        # DIA results do not have spec_idx/scan_num in psm_df, nothing to merge
        return psm_df

    def _prepare_matching_dfs(self):
        fragment_mz_df = self.get_fragment_mz_df()
        fragment_mz_df = pd.concat(
            [fragment_mz_df] * self.max_spec_per_query, ignore_index=True
        )
        if self.use_ppm:
            self.all_frag_mz_tols = fragment_mz_df.values * self.tolerance * 1e-6
        else:
            self.all_frag_mz_tols = np.full_like(fragment_mz_df.values, self.tolerance)

        psm_df_list = []
        len_frags = len(fragment_mz_df) // self.max_spec_per_query
        for i in range(self.max_spec_per_query):
            psm_df = self.psm_df.copy()
            psm_df["frag_start_idx"] = psm_df.frag_start_idx + i * len_frags
            psm_df["frag_stop_idx"] = psm_df.frag_stop_idx + i * len_frags
            psm_df_list.append(psm_df)
        self.psm_df = pd.concat(psm_df_list, ignore_index=True)

        matched_intensity_df = pd.DataFrame(
            np.zeros_like(fragment_mz_df.values, dtype=PEAK_INTENSITY_DTYPE),
            columns=fragment_mz_df.columns,
        )

        matched_mz_err_df = pd.DataFrame(
            np.zeros_like(fragment_mz_df.values, dtype=PEAK_MZ_DTYPE),
            columns=fragment_mz_df.columns,
        )

        return (fragment_mz_df, matched_intensity_df, matched_mz_err_df)

    def _match_ms2_one_raw_numba(
        self, raw_name: str, psm_df_one_raw: pd.DataFrame
    ) -> pd.DataFrame:
        """
        Internal method to extract peak information with numba as backend.

        Parameters
        ----------
        raw_name : str
            The raw name of the raw file. `psm_df_one_raw` dataframe should also
            contain the same raw name in `raw_name` column.
        psm_df_one_raw : pd.DataFrame
            The dataframe for PSMs.

        Returns
        -------
        pd.DataFrame
            `psm_df_one_raw`
        """
        psm_df_one_raw = psm_df_one_raw.reset_index(drop=True)

        if raw_name in self._ms_file_dict:
            raw_data = load_ms_data(
                self._ms_file_dict[raw_name],
                self._ms_file_type,
                process_count=self.ms_loader_thread_num,
            )
            if raw_data is None:
                return

            psm_origin_len = len(psm_df_one_raw) // self.max_spec_per_query

            grouper = NormalDIAGrouper(raw_data)

            psm_groups = grouper.assign_dia_groups(
                psm_df_one_raw.precursor_mz.values[:psm_origin_len]
            )

            all_spec_idxes = np.full(len(psm_df_one_raw), -1, dtype=np.int32)

            for dia_group, group_df in grouper.dia_group_dfs:
                psm_idxes = psm_groups[dia_group]
                if len(psm_idxes) == 0:
                    continue
                psm_idxes = np.array(psm_idxes, dtype=np.int32)
                spec_idxes = find_dia_spec_idxes_same_window(
                    group_df.rt.values,
                    psm_df_one_raw.rt.values[psm_idxes],
                    max_spec_per_query=self.max_spec_per_query,
                )
                for i in range(spec_idxes.shape[-1]):
                    all_spec_idxes[psm_idxes + psm_origin_len * i] = spec_idxes[:, i]

                match_one_raw_with_numba(
                    all_spec_idxes,
                    psm_df_one_raw.frag_start_idx.values,
                    psm_df_one_raw.frag_stop_idx.values,
                    self.fragment_mz_df.values,
                    self.all_frag_mz_tols,
                    raw_data.peak_df.mz.values,
                    raw_data.peak_df.intensity.values,
                    group_df.peak_start_idx.values,
                    group_df.peak_stop_idx.values,
                    self.matched_intensity_df.values,
                    self.matched_mz_err_df.values,
                    self.match_closest,
                )
        else:
            print(f"`{raw_name}` is not found in ms_file_dict.")
            return
        return psm_df_one_raw

    def match_ms2_multi_raw(
        self,
        psm_df: pd.DataFrame,
        ms_files: Tuple[dict, list],
        ms_file_type: str = "alpharaw_hdf",
        process_num: int = 8,
    ) -> Tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
        """
        Match peaks for the given `psm_df` against the corresponding MS spectrum files.

        Parameters
        ----------
        psm_df : pd.DataFrame
            Peptide-spectrum matches in alphabase dataframe format.
        ms_files : Tuple[dict, list]
            The absolute or relative paths of MS files.
            if the type is `dict`, the format will be
            `{'raw_name1': 'raw_name1.raw', ...}` if `ms_file_type` is `thermo_raw`.
        ms_file_type : str, optional
            MS file type that is already registered in
            :obj:`alpharaw.ms_data_base.ms_reader_provider`.
            By default "alpharaw_hdf".
        process_num : int, optional
            Match peaks by using multiprocessing, by default 8

        Returns
        -------
        Tuple
            pd.DataFrame: the `psm_df`.

            pd.DataFrame: fragment m/z dataframe in alphabase format.

            pd.DataFrame: the matched fragment intensity dataframe in alphabase format.

            pd.DataFrame: the matched mass error in the same dataframe format.
        """
        if isinstance(ms_files, list):
            ms_files = parse_ms_files_to_dict(ms_files)
        psm_df = psm_df[psm_df.raw_name.isin(ms_files)].reset_index(drop=True)
        super().match_ms2_multi_raw(psm_df, ms_files, ms_file_type, process_num)

        return (
            self.psm_df,
            self.fragment_mz_df,
            self.matched_intensity_df,
            self.matched_mz_err_df,
        )
