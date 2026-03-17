"""DIA grouping for normal DIA data, which is not acquired with ion mobility separation. Moved from alpharaw."""
import typing
from collections import defaultdict

import numpy as np
from alphatims.alphatims_wrapper import convert_to_alphatims, remove_unused_peaks
from alphatims.bruker import TimsTOF

from alpharaw.ms_data_base import MSData_Base


class NormalDIAGrouper:
    def __init__(self, ms_data: MSData_Base):
        self.ms_data = ms_data
        self.ms_data.spectrum_df["dia_group"] = ms_data.spectrum_df.precursor_mz.astype(
            int
        )

        self.dia_group_dfs = self.ms_data.spectrum_df.groupby("dia_group")
        self.dia_isolation_dict = {}
        for dia_group, df in self.dia_group_dfs:
            if dia_group == -1:
                continue
            self.dia_isolation_dict[dia_group] = (
                df.isolation_lower_mz.values[0],
                df.isolation_upper_mz.values[0],
            )
        self.dia_groups = np.sort(list(self.dia_isolation_dict.keys()))

    def get_ms_data_for_a_group(
        self,
        dia_group: int = -1,
        return_alpharaw_data: bool = True,
        return_alphatims_data: bool = True,
    ) -> typing.Union[MSData_Base, TimsTOF, typing.Tuple[MSData_Base, TimsTOF]]:
        """Get compressed MS data for isolation window `dia_group`.

        Args:
            dia_group (int, optional): The DIA group, -1 means ms1. Defaults to -1.
            return_alphatims_data (bool, optional): If return `MSData_Base`. Defaults to True
            return_alphatims_data (bool, optional): If return alphatims object. Defaults to True.

        Returns:
            MSData_Base: Compressed MS data, if `return_alpharaw_data==True`
            TimsTOF: Alphatims object for the window, if `return_alphatims_data==True`
        """

        spec_df = self.dia_group_dfs.get_group(dia_group)

        if return_alphatims_data:
            ms_data, ms_tims = convert_to_alphatims(
                spec_df, self.ms_data.peak_df, dda=False
            )
            if return_alpharaw_data:
                return ms_data, ms_tims
            else:
                return ms_tims
        else:
            ms_data = MSData_Base()

            spec_df, peak_df = remove_unused_peaks(spec_df, self.ms_data.peak_df)

            ms_data.spectrum_df = spec_df
            ms_data.peak_df = peak_df
            return ms_data

    def assign_dia_groups(
        self, precursor_mzs
    ) -> typing.DefaultDict[typing.List, typing.List]:
        dia_precursor_groups = defaultdict(list)
        for i, mz in enumerate(precursor_mzs):
            i_group = np.searchsorted(self.dia_groups, int(mz))
            if i_group == 0:
                i_group = 1
            elif i_group == len(self.dia_groups):
                i_group -= 1
            dia_group = self.dia_groups[i_group]
            if dia_group in self.dia_isolation_dict:
                isolation_lower, isolation_upper = self.dia_isolation_dict[dia_group]
                if mz >= isolation_lower and mz <= isolation_upper:
                    dia_precursor_groups[dia_group].append(i)
                    continue
            dia_group = self.dia_groups[i_group - 1]
            if dia_group in self.dia_isolation_dict:
                isolation_lower, isolation_upper = self.dia_isolation_dict[dia_group]
                if mz >= isolation_lower and mz <= isolation_upper:
                    dia_precursor_groups[dia_group].append(i)
        return dia_precursor_groups
