import os
import psutil

import pandas as pd
import numpy as np
from typing import Union

from alphabase.peptide.fragment import get_charged_frag_types

from peptdeep.settings import global_settings
from peptdeep.protein.fasta import PredictSpecLibFasta
from peptdeep.spec_lib.translate import (
    speclib_to_single_df, mod_to_unimod_dict,
    translate_to_tsv
)

from peptdeep.pretrained_models import ModelManager
from peptdeep.utils import logging

class PredictLibraryMakerBase(object):
    """
    Base class to predict libraries
    """
    def __init__(self, 
        model_manager:ModelManager = None,
    ):
        lib_settings = global_settings['library']
        self.spec_lib = PredictSpecLibFasta(
            model_manager=model_manager,
            charged_frag_types = get_charged_frag_types(
                lib_settings['frag_types'],
                lib_settings['max_frag_charge'],
            ),
            protease = lib_settings['fasta']['protease'],
            max_missed_cleavages = lib_settings['fasta']['max_miss_cleave'],
            peptide_length_min = lib_settings['min_peptide_len'],
            peptide_length_max = lib_settings['max_peptide_len'],
            precursor_charge_min = lib_settings['min_precursor_charge'],
            precursor_charge_max = lib_settings['max_precursor_charge'],
            precursor_mz_min = lib_settings['min_precursor_mz'], 
            precursor_mz_max = lib_settings['max_precursor_mz'],
            var_mods = lib_settings['var_mods'],
            min_var_mod_num = lib_settings['min_var_mod_num'],
            max_var_mod_num = lib_settings['max_var_mod_num'],
            fix_mods = lib_settings['fix_mods'],
            labeling_channels = lib_settings['labeling_channels'],
            special_mods = lib_settings['special_mods'],
            min_special_mod_num = lib_settings['min_special_mod_num'],
            max_special_mod_num = lib_settings['max_special_mod_num'],
            special_mods_cannot_modify_pep_n_term = lib_settings['special_mods_cannot_modify_pep_n_term'],
            special_mods_cannot_modify_pep_c_term = lib_settings['special_mods_cannot_modify_pep_c_term'],
            decoy = lib_settings['decoy'],
            include_contaminants=lib_settings['fasta']['add_contaminants'],
            I_to_L=False,
            generate_precursor_isotope=lib_settings['generate_precursor_isotope'],
            rt_to_irt=lib_settings['rt_to_irt'],
    )

    def _check_df(self)->str:
        pass

    def _input(self, _input):
        """Virtual method to be re-implemented by sub-classes"""
        raise NotImplementedError("All sub-classes must re-implement '_input()' method")

    def _predict(self):
        self.spec_lib.predict_all()

    @property
    def precursor_df(self)->pd.DataFrame:
        return self.spec_lib.precursor_df

    @property
    def fragment_intensity_df(self)->pd.DataFrame:
        return self.spec_lib.fragment_intensity_df

    @property
    def fragment_mz_df(self)->pd.DataFrame:
        return self.spec_lib.fragment_mz_df

    def make_library(self, _input):
        """Predict a library for the `_input`, 
        this function runs the following methods.

        - self._input(_input)
        - self._check_df()
        - self._predict()

        Parameters
        ----------
        _input
            _input file or source

        Raises
        ------
        ValueError
            ValueError for some reasons
        """
        logging.info("Generating the spectral library ...")
        try:
            self._input(_input)
            self._check_df()
            self._predict()

            logging.info(
                'Predicting the spectral library with '
                f'{len(self.precursor_df)} precursors '
                f'and {np.prod(self.fragment_mz_df.values.shape, dtype=float)*(1e-6):.2f}M fragments '
                f'used {psutil.Process(os.getpid()).memory_info().rss/1024**3:.4f} GB memory'
            )
        except ValueError as e:
            raise e
    
    def translate_to_tsv(self, 
        tsv_path:str, 
        translate_mod_dict:dict=None
    ):
        """Translate the predicted DataFrames into a TSV file
        """
        logging.info(f"Translating to {tsv_path} for DiaNN/Spectronaut...")
        lib_settings = global_settings['library']

        if 'proteins' not in self.spec_lib._precursor_df.columns:
            self.spec_lib.append_protein_name()
        
        translate_to_tsv(
            self.spec_lib, 
            tsv_path,
            keep_k_highest_fragments=lib_settings['output_tsv'][
                'keep_higest_k_peaks'
            ],
            min_frag_intensity=lib_settings['output_tsv'][
                'min_relative_intensity'
            ],
            min_frag_mz=lib_settings['output_tsv'][
                'min_fragment_mz'
            ],
            max_frag_mz=lib_settings['output_tsv'][
                'max_fragment_mz'
            ],
            batch_size=lib_settings['output_tsv'][
                'translate_batch_size'
            ],
            translate_mod_dict=translate_mod_dict,
        )
    
    def translate_library(self, 
        translate_mod_dict:dict=None
    )->pd.DataFrame:
        """Translate predicted DataFrames into 
        a single DataFrame in SWATH library format
        """
        logging.info("Translating library for DiaNN/Spectronaut...")
        lib_settings = global_settings['library']

        if 'proteins' not in self.spec_lib._precursor_df.columns:
            self.spec_lib.append_protein_name()
        
        return speclib_to_single_df(
            self.spec_lib, 
            translate_mod_dict=translate_mod_dict,
            keep_k_highest_fragments=lib_settings['output_tsv'][
                'keep_higest_k_peaks'
            ],
            min_frag_intensity=lib_settings['output_tsv'][
                'min_relative_intensity'
            ],
            min_frag_mz=lib_settings['output_tsv'][
                'min_fragment_mz'
            ],
            max_frag_mz=lib_settings['output_tsv'][
                'max_fragment_mz'
            ],
        )

class PrecursorLibraryMaker(PredictLibraryMakerBase):
    """For input dataframe of charged modified sequences"""
    def _input(self, precursor_df:pd.DataFrame):
        self.spec_lib._precursor_df = precursor_df
        self.spec_lib.add_peptide_labeling()
        self.spec_lib.append_decoy_sequence()
    
    def _check_df(self):
        if 'charge' not in self.spec_lib.precursor_df.columns:
            raise ValueError('self.spec_lib.precursor_df must contain the "charge" column.')
        (
            self.spec_lib.precursor_df['charge']
        ) = self.spec_lib.precursor_df['charge'].astype(np.int8)

        if (
            'mods' not in self.spec_lib.precursor_df.columns or 
            'mod_sites' not in self.spec_lib.precursor_df.columns
        ):
            self.spec_lib.precursor_df['mods'] = ''
            self.spec_lib.precursor_df['mod_sites'] = ''
        else:
            (
                self.spec_lib.precursor_df['mods']
            ) = self.spec_lib.precursor_df['mods'].astype(str)
            (
                self.spec_lib.precursor_df['mod_sites']
            ) = self.spec_lib.precursor_df['mod_sites'].astype(str)

        self.spec_lib.protein_df = pd.DataFrame()

class PeptideLibraryMaker(PrecursorLibraryMaker):
    """For input dataframe of modified sequences"""
    def _input(self, peptide_df:pd.DataFrame):
        self.spec_lib._precursor_df = peptide_df
        self.spec_lib.append_decoy_sequence()
        self.spec_lib.add_peptide_labeling()
        self.spec_lib.add_charge()

class SequenceLibraryMaker(PeptideLibraryMaker):
    """For input dataframe of AA sequences"""
    def _input(self, sequence_df:pd.DataFrame):
        self.spec_lib._precursor_df = sequence_df
        self.spec_lib.append_decoy_sequence()
        self.spec_lib.add_modifications()
        self.spec_lib.add_special_modifications()
        self.spec_lib.add_peptide_labeling()
        self.spec_lib.add_charge()

class FastaLibraryMaker(PredictLibraryMakerBase):
    """For fasta or a list of fasta files"""
    def _input(self, fasta:Union[str,list]):
        self.spec_lib.get_peptides_from_fasta(fasta)
        self.spec_lib.append_decoy_sequence()
        self.spec_lib.add_modifications()
        self.spec_lib.add_special_modifications()
        self.spec_lib.add_peptide_labeling()
        self.spec_lib.add_charge()

class LibraryMakerProvider:
    """
    Factory class for library makers
    """
    def __init__(self):
        self.library_maker_dict = {}

    def register_maker(self, maker_name:str, maker_class):
        self.library_maker_dict[maker_name.lower()] = maker_class

    def get_maker(self, maker_name:str, *,
        model_manager = None,
    )->PredictLibraryMakerBase:
        maker_name = maker_name.lower()
        if maker_name in self.library_maker_dict:
            return self.library_maker_dict[maker_name](model_manager)
        else:
            raise ValueError(f'library maker "{maker_name}" is not registered.')

library_maker_provider = LibraryMakerProvider()
library_maker_provider.register_maker('precursor_table', PrecursorLibraryMaker)
library_maker_provider.register_maker('precursor_library', PrecursorLibraryMaker)
library_maker_provider.register_maker('peptide_table', PeptideLibraryMaker)
library_maker_provider.register_maker('peptide_library', PeptideLibraryMaker)
library_maker_provider.register_maker('sequence_table', SequenceLibraryMaker)
library_maker_provider.register_maker('sequence_library', SequenceLibraryMaker)
library_maker_provider.register_maker('fasta', FastaLibraryMaker)
library_maker_provider.register_maker('fasta_library', FastaLibraryMaker)
