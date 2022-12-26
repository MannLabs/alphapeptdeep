import pandas as pd
import numpy as np
import torch

from alphabase.peptide.precursor import (
    calc_precursor_isotope_mp, calc_precursor_isotope
)
from alphabase.spectral_library.base import SpecLibBase
from alphabase.spectral_library.flat import SpecLibFlat
from alphabase.peptide.fragment import (
    flatten_fragments, concat_precursor_fragment_dataframes
)

from peptdeep.pretrained_models import ModelManager
from peptdeep.settings import global_settings
from peptdeep.utils import logging
from peptdeep.utils import process_bar

lib_settings = global_settings['library']
model_mgr_settings = global_settings['model_mgr']

class PredictSpecLib(SpecLibBase):
    """PredictSpecLib

    Parameters
    ----------
    model_manager : ModelManager, optional
        `ModelManager`, by default None

    charged_frag_types : list, optional
        Charged fragment types, by default ['b_z1','b_z2','y_z1','y_z2']

    precursor_mz_min : float, optional
        precursor_mz_min, by default 400.0

    precursor_mz_max : float, optional
        precursor_mz_max, by default 2000.0

    generate_precursor_isotope : bool, optional
        If calculate isotope masses and relative intensities for precursors

    decoy : str, optional
        Decoy choice, see `alphabase.spec_lib.decoy_library`, 
        by default 'pseudo_reverse'
    """
    def __init__(self,
        model_manager: ModelManager = None,
        charged_frag_types = ['b_z1','b_z2','y_z1','y_z2'],
        precursor_mz_min:float = 400.0, 
        precursor_mz_max:float = 2000.0,
        generate_precursor_isotope:bool = False,
        decoy:str = 'pseudo_reverse'
    ):
        super().__init__(
            charged_frag_types,
            precursor_mz_min=precursor_mz_min,
            precursor_mz_max=precursor_mz_max,
            decoy = decoy
        )
        self.generate_precursor_isotope = generate_precursor_isotope
        self.verbose = True
        if model_manager is None:
            self.model_manager = ModelManager(
                mask_modloss=False
            )
        else:
            self.model_manager = model_manager

        self._precursor_df = pd.DataFrame()
        self._fragment_intensity_df = pd.DataFrame()
        self._fragment_mz_df = pd.DataFrame()

        self.mp_predict_batch_size:int = 100000
        self.use_multiprocessing:bool = model_mgr_settings['predict']['multiprocessing']
        self.mp_predict_process_num:int = global_settings['thread_num']

    def set_precursor_and_fragment(self,
        *,
        precursor_df: pd.DataFrame,
        fragment_mz_df: pd.DataFrame,
        fragment_intensity_df: pd.DataFrame,
    ):
        self._precursor_df = precursor_df
        self._fragment_intensity_df = fragment_intensity_df
        self._fragment_mz_df = fragment_mz_df

        self._fragment_mz_df.drop(columns=[
            col for col in self._fragment_mz_df.columns 
            if col not in self.charged_frag_types
        ], inplace=True)

        self._fragment_intensity_df.drop(columns=[
            col for col in self._fragment_intensity_df.columns 
            if col not in self.charged_frag_types
        ], inplace=True)

    def rt_to_irt_pred(self):
        """ Add 'irt_pred' into columns based on 'rt_pred' """
        return self.model_manager.rt_model.add_irt_column_to_precursor_df(self._precursor_df)

    def predict_all(self, 
        min_required_precursor_num_for_mp:int=2000,
    ):
        """
        1. Predict RT/IM/MS2 for self._precursor_df
        2. Calculate isotope information in self._precursor_df
        """
        self.calc_precursor_mz()
        if self.generate_precursor_isotope:
            if self.verbose:
                logging.info('Calculating precursor isotope distributions ...')
            if len(self.precursor_df) < min_required_precursor_num_for_mp:
                self._precursor_df = calc_precursor_isotope(
                    self._precursor_df
                )
            else:
                self._precursor_df = calc_precursor_isotope_mp(
                    self._precursor_df, process_bar=process_bar
                )
        if self.verbose:
            logging.info('Predicting RT/IM/MS2 ...')
        res = self.model_manager.predict_all(
            self._precursor_df,
            predict_items=['rt','mobility','ms2'],
            frag_types=self.charged_frag_types,
            min_required_precursor_num_for_mp=min_required_precursor_num_for_mp,
            multiprocessing=self.use_multiprocessing,
            mp_batch_size=self.mp_predict_batch_size,
            process_num=self.mp_predict_process_num,
        )
        self.set_precursor_and_fragment(**res)
        if self.verbose:
            logging.info('End Predicting RT/IM/MS2')
        

class PredictSpecLibFlat(SpecLibFlat):
    """ 
    Flatten the predicted spectral library, the key feature is to 
    predict and flatten fragments in batch with `predict_and_parse_lib_in_batch()`

    Parameters
    ----------
    min_fragment_intensity : float, optional
        minimal intensity to keep, by default 0.001
    keep_top_k_fragments : int, optional
        top k highest peaks to keep, by default 1000
    """
    def __init__(self, 
        min_fragment_intensity:float = 0.001,
        keep_top_k_fragments:int = 1000,
        custom_fragment_df_columns:list = [
            'type','number','position','charge','loss_type'
        ],
        **kwargs,
    ):
        super().__init__(
            min_fragment_intensity=min_fragment_intensity,
            keep_top_k_fragments=keep_top_k_fragments,
            custom_fragment_df_columns=custom_fragment_df_columns
        )

    def predict_and_parse_lib_in_batch(self, 
        predict_lib:PredictSpecLib, 
        batch_size:int = 200000
    ):
        """Predict and flatten fragments in batch

        Parameters
        ----------
        predict_lib : PredictSpecLib
            spectral library to be predicted and flatten
        batch_size : int, optional
            the batch size, by default 200000
        """
        if len(predict_lib.precursor_df) <= batch_size:
            predict_lib.predict_all()
            self.parse_base_library(predict_lib)
        else:
            predict_lib.verbose = False
            predict_lib.refine_df()
            precursor_df = predict_lib.precursor_df
            precursor_df_list = []
            fragment_df_list = []
            for i in range(0, len(precursor_df), batch_size):
                predict_lib._precursor_df = precursor_df.iloc[i:i+batch_size].copy()
                predict_lib.predict_all()
                df, frag_df = flatten_fragments(
                    predict_lib.precursor_df,
                    predict_lib.fragment_mz_df,
                    predict_lib.fragment_intensity_df,
                    min_fragment_intensity = self.min_fragment_intensity,
                    keep_top_k_fragments = self.keep_top_k_fragments,
                    custom_columns=self.custom_fragment_df_columns
                )
                precursor_df_list.append(df)
                fragment_df_list.append(frag_df)
            predict_lib._precursor_df = precursor_df
            self._precursor_df, self._fragment_df = concat_precursor_fragment_dataframes(
                precursor_df_list, fragment_df_list
            )


