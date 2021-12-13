# AUTOGENERATED! DO NOT EDIT! File to edit: nbdev_nbs/spec_lib/predict_lib.ipynb (unless otherwise specified).

__all__ = ['PredictLib']

# Cell
import pandas as pd

from alphabase.spectrum_library.library_base import SpecLibBase
from alphabase.peptide.fragment import update_precursor_mz
from alphadeep.pretrained_models import AlphaDeepModels

class PredictLib(SpecLibBase):
    def __init__(self,
        models: AlphaDeepModels,
        charged_frag_types, #['b_z1','b_z2','y_z1','y_z2', ...]
        min_frag_mz = 50, max_frag_mz = 2000,
        min_precursor_mz = 400, max_precursor_mz = 2000,
    ):
        super().__init__(
            charged_frag_types,
            min_frag_mz=min_frag_mz,
            max_frag_mz=max_frag_mz,
            min_precursor_mz=min_precursor_mz,
            max_precursor_mz=max_precursor_mz
        )
        self.models = models

        self.intensity_factor = 1
        self.verbose = True

        self._precursor_df = pd.DataFrame()
        self._fragment_intensity_df = pd.DataFrame()
        self._fragment_mz_df = pd.DataFrame()

    @property
    def precursor_df(self):
        return self._precursor_df

    @precursor_df.setter
    def precursor_df(self, df):
        self._precursor_df = df
        self._init_precursor_df()

    def _init_precursor_df(self):
        self._precursor_df['nAA'] = self._precursor_df['sequence'].str.len()
        self._precursor_df['mod_sites'] = self._precursor_df['mod_sites'].astype('U')
        self._precursor_df['charge'] = self._precursor_df['charge'].astype(int)
        if 'precursor_mz' not in self._precursor_df.columns:
            update_precursor_mz(self._precursor_df)

    def predict_rt_ccs(self):
        # add 'rt_pred' and 'irt_pred' into columns
        self._precursor_df = self.models.rt_model.predict(
            self._precursor_df, verbose=self.verbose
        )
        self.models.rt_model.rt_to_irt_pred(self._precursor_df)
        # add 'ccs_pred' and 'mobility_pred' into columns
        self._precursor_df = self.models.ccs_model.predict(
            self._precursor_df, verbose=self.verbose
        )
        self.models.ccs_model.ccs_to_mobility_pred(
            self._precursor_df
        )

    def load_fragment_intensity_df(self, **kwargs):
        if len(self._fragment_mz_df) == 0:
            self.calc_fragment_mz_df()

        frag_inten_df = self.models.ms2_model.predict(
            self._precursor_df,
            reference_frag_df=self._fragment_mz_df,
            verbose=self.verbose,
        )

        charged_frag_list = []
        for frag_type in self._fragment_mz_df.columns.values:
            if frag_type in frag_inten_df:
                charged_frag_list.append(frag_type)
        self._fragment_mz_df = self._fragment_mz_df[
            charged_frag_list
        ]
        self._fragment_intensity_df = frag_inten_df[
            charged_frag_list
        ]*self.intensity_factor
        self._fragment_intensity_df[self._fragment_mz_df==0] = 0

