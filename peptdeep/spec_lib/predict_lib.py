import pandas as pd
import numpy as np
import torch
import tqdm

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

model_mgr_settings = global_settings['model_mgr']

class PredictSpecLib(SpecLibBase):
    def __init__(self,
        model_manager: ModelManager = None,
        charged_frag_types = ['b_z1','b_z2','y_z1','y_z2'],
        precursor_mz_min:float = 400.0, 
        precursor_mz_max:float = 2000.0,
        decoy:str = 'pseudo_reverse',
        rt_to_irt:bool = False,
        generate_precursor_isotope:bool = False,
    ):
        """
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

        decoy : str, optional
            Decoy choice, see `alphabase.spec_lib.decoy_library`, 
            by default 'pseudo_reverse'

        rt_to_irt : bool, optional
            Convert predicted RT to iRT values, by default False

        generate_precursor_isotope : bool, optional
            Generate precursor isotopes, defaults to False
        """
        super().__init__(
            charged_frag_types,
            precursor_mz_min=precursor_mz_min,
            precursor_mz_max=precursor_mz_max,
            decoy = decoy
        )
        if model_manager is None:
            self.model_manager = ModelManager(
                mask_modloss=True
            )
        else:
            self.model_manager = model_manager

        self._precursor_df = pd.DataFrame()
        self._fragment_intensity_df = pd.DataFrame()
        self._fragment_mz_df = pd.DataFrame()

        self.mp_predict_batch_size:int = 100000
        self.rt_to_irt = rt_to_irt
        self.generate_precursor_isotope = generate_precursor_isotope

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

    def translate_rt_to_irt_pred(self):
        """ Add 'irt_pred' into columns based on 'rt_pred' """
        return self.model_manager.rt_model.add_irt_column_to_precursor_df(self._precursor_df)

    def predict_all(self, 
        min_required_precursor_num_for_mp:int=2000,
        predict_items:list = ['rt','mobility','ms2'],
    ):
        """
        1. Predict RT/IM/MS2 for self._precursor_df
        2. Calculate isotope information in self._precursor_df
        """
        self.calc_precursor_mz()
        if self.generate_precursor_isotope:
            if self.model_manager.verbose:
                logging.info('Calculating precursor isotope distributions ...')
            if len(self.precursor_df) < min_required_precursor_num_for_mp:
                self._precursor_df = calc_precursor_isotope(
                    self._precursor_df
                )
            else:
                self._precursor_df = calc_precursor_isotope_mp(
                    self._precursor_df, process_bar=process_bar
                )
        if self.model_manager.verbose:
            logging.info(f'Predicting RT/IM/MS2 for {len(self._precursor_df)} precursors ...')
        res = self.model_manager.predict_all(
            self._precursor_df,
            predict_items=predict_items,
            frag_types=self.charged_frag_types,
            min_required_precursor_num_for_mp=min_required_precursor_num_for_mp,
            multiprocessing=model_mgr_settings['predict']['multiprocessing'],
            mp_batch_size=self.mp_predict_batch_size,
            process_num=global_settings['thread_num'],
        )
        self.set_precursor_and_fragment(**res)
        if self.rt_to_irt and 'rt_pred' in self._precursor_df.columns:
            self.translate_rt_to_irt_pred()
        if self.model_manager.verbose:
            logging.info('End predicting RT/IM/MS2')
        

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
        logging.info(f"Flattening {len(predict_lib.precursor_df)} precursors in batch size {batch_size} ...")
        if len(predict_lib.precursor_df) <= batch_size:
            predict_lib.predict_all()
            self.parse_base_library(predict_lib)
        else:
            predict_lib.model_manager.verbose = False
            predict_lib.refine_df()
            df = predict_lib.precursor_df
            precursor_df_list = []
            fragment_df_list = []
            for i in tqdm.tqdm(range(0, len(df), batch_size)):
                predict_lib._precursor_df = df.iloc[i:i+batch_size].copy()
                predict_lib.predict_all()
                flat_df, frag_df = flatten_fragments(
                    predict_lib.precursor_df,
                    predict_lib.fragment_mz_df,
                    predict_lib.fragment_intensity_df,
                    min_fragment_intensity = self.min_fragment_intensity,
                    keep_top_k_fragments = self.keep_top_k_fragments,
                    custom_columns=self.custom_fragment_df_columns
                )
                precursor_df_list.append(flat_df)
                fragment_df_list.append(frag_df)
            predict_lib._precursor_df = df
            self._precursor_df, self._fragment_df = concat_precursor_fragment_dataframes(
                precursor_df_list, fragment_df_list
            )


