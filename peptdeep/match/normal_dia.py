"""DIA grouping for normal DIA data, which is not acquired with ion mobility separation. Moved from alpharaw."""

import typing
from collections import defaultdict

import numpy as np

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
