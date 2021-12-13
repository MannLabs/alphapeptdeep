# AUTOGENERATED! DO NOT EDIT! File to edit: nbdev_nbs/psm_reader/psm_reader.ipynb (unless otherwise specified).

__all__ = ['translate_other_modification', 'keep_modifications', 'PSMReaderBase', 'PSMReader_w_FragBase',
           'PSMReaderProvider', 'psm_reader_provider', 'PSMwFragReaderProvider', 'psm_w_frag_reader_provider']

# Cell
import typing
import numpy as np
import pandas as pd

from alphabase.peptide.fragment import get_charged_frag_types

def translate_other_modification(
    mod_str: str,
    mod_dict: dict
)->str:
    '''
    Translate modifications in `mod_str` to the other
    format mapped by mod_dict.
    Args:
        mod_str (str): mod list in str format, seperated by ';',
            e.g. ModA;ModB
        mod_dict (dict): translate mod dict from others to AlphaBase,
            e.g. for pFind, key='Phospho[S]', value='Phospho@S'
    Returns:
        str: new mod list in str format seperated by ';' if all
             modifications are in `mod_dict` else pd.NA.
    '''
    if not mod_str: return ""
    ret_mods = []
    for mod in mod_str.split(';'):
        if mod in mod_dict:
            ret_mods.append(mod_dict[mod])
        else:
            return pd.NA
    return ";".join(ret_mods)

def keep_modifications(
    mod_str: str,
    mod_set: set
)->str:
    '''
    Check if modifications in `mod_str` in `mod_set`
    Args:
        mod_str (str): mod list in str format, seperated by ';',
            e.g. Oxidation@M;Phospho@S.
        mod_set (set): mod set to check
    Returns:
        str: `mod_str` if all modifications are in mod_set
             else pd.NA.
    '''
    if not mod_str: return ""
    for mod in mod_str.split(';'):
        if not mod in mod_set:
            return pd.NA
    return mod_str


class PSMReaderBase(object):
    def __init__(self,
    ):
        # modification_convert_dict=dict[str, str]:
        #     key:   mod names of other search engines
        #     value: mod names in AlphaBase
        # It is used to convert mods of other engines
        # to AlphaBase format. Different search engines
        # have different mod names.

        self.modification_convert_dict = {}

        self.column_mapping = {
            'sequence': 'NakedSequence',
            # AlphaBase does not need 'modified_sequence',
            # but it will get 'mods', 'mod_sites' from it
            'modified_sequence': 'ModifiedSequence',
            'charge': 'Charge',
            # If the value is a list, check if one of the columns exist
            # and get 'proteins' from that column
            'proteins':['Proteins','UniprotIDs'],
            'uniprot_ids':'UniprotIds',
            # Similar to 'proteins'
            'genes': ['Genes','Gene Names','Gene names'],
        } # Add more columns for sub-classes of different tasks

        self._psm_df:pd.DataFrame = None
        self.keep_all_psm = False

    @property
    def psm_df(self):
        return self._psm_df

    def _load_file(self, filename:str)->pd.DataFrame:
        """
        Load original dataframe from PSM filename.
        Different search engines may store PSMs in different ways:
        tsv, csv, HDF, XML, ...

        Args:
            filename (str): psm filename

        Raises:
            NotImplementedError: Sub-classes must re-implement this method

        Returns:
            pd.DataFrame: dataframe loaded
        """
        raise NotImplementedError(
            f'"{self.__class__}" must re-implement "_load_file()"'
        )

    def _translate_columns(self, origin_df:pd.DataFrame):
        """
        Translate the dataframe from other search engines
        to AlphaBase format

        Args:
            origin_df (pd.DataFrame): df of other search engines
        """
        self._psm_df = pd.DataFrame()
        for col, map_col in self.column_mapping.items():
            if isinstance(map_col, str):
                if map_col in origin_df.columns:
                    self._psm_df[col] = origin_df[map_col]
                # else:
                #     self._psm_df[col] = pd.NA
            else:
                for other_col in map_col:
                    if other_col in origin_df.columns:
                        self._psm_df[col] = origin_df[other_col]
                        break
                # if col not in self._psm_df.columns:
                #     self._psm_df[col] = pd.NA

    def _translate_modifications(self):
        '''
        Translate modifications to AlphaBase format.

        Raises: KeyError if `mod` in `mod_names` is
            not in `self.modification_convert_dict`
        '''
        self._psm_df.mods = self._psm_df.mods.apply(
            translate_other_modification,
            mod_dict=self.modification_convert_dict
        )

    def _post_process(self,
        filename:str, origin_df:pd.DataFrame
    ):
        """
        Post processing after loading and translate everything.
        Here, we remove unknown modifications and perform other post processings,
        e.g. loading fragments for AlphaQuant or AlphaDeep

        Args:
            filename (str): psm filename
            origin_df (pd.DataFrame): the loaded original df
        """
        self._psm_df['nAA'] = self._psm_df.sequence.str.len()
        origin_df = origin_df[
            ~self._psm_df['mods'].isna()
        ].reset_index(drop=True)

        self._psm_df = self._psm_df[
            ~self._psm_df['mods'].isna()
        ].reset_index(drop=True)

    def load(self, filename):
        origin_df = self._load_file(filename)
        self._translate_columns(origin_df)
        self._translate_modifications()
        self._post_process(filename, origin_df)

    def normalize_rt_by_raw_name(self):
        if (
            not 'raw_name' in self.psm_df.columns
            or not 'rt_norm' in self.psm_df.columns
        ):
            return
        for raw_name, df_group in self.psm_df.groupby('raw_name'):
            self.psm_df.loc[
                df_group.index,'rt_norm'
            ] = df_group.rt_norm / df_group.rt_norm.max()

    def filter_psm_by_modifications(self, include_mod_list = [
        'Oxidation@M','Phospho@S','Phospho@T','Phospho@Y','Acetyl@Protein N-term'
    ]):
        '''
            Only keeps peptides with modifications in `include_mod_list`.
        '''
        mod_set = set(include_mod_list)
        self._psm_df.mods = self._psm_df.mods.apply(keep_modifications, mod_set=mod_set)

        self._psm_df.dropna(
            subset=['mods'], inplace=True
        )
        self._psm_df.reset_index(drop=True, inplace=True)

class PSMReader_w_FragBase(PSMReaderBase):
    '''
    Read PSMs and fragments
    '''
    def __init__(self,
        frag_types=['b','y','b_modloss','y_modloss'],
        max_frag_charge=2,
        frag_tol=20, frag_ppm=True,
    ):
        super().__init__()

        self.charged_frag_types = get_charged_frag_types(
            frag_types, max_frag_charge
        )
        self._fragment_intensity_df:pd.DataFrame = pd.DataFrame(
            columns=self.charged_frag_types
        )

        self.frag_tol = frag_tol
        self.frag_ppm = frag_ppm

    @property
    def fragment_intensity_df(self):
        return self._fragment_intensity_df


# Cell
class PSMReaderProvider:
    def __init__(self):
        self.reader_dict = {}

    def register_reader(self, reader_name, reader_class):
        # for example, we can register the MSFragger reader
        self.reader_dict[reader_name.lower()] = reader_class

    def get_reader(self, reader_name,
    )->PSMReaderBase:
        return self.reader_dict[reader_name.lower()]()

psm_reader_provider = PSMReaderProvider()

# Cell
class PSMwFragReaderProvider:
    def __init__(self):
        self.reader_dict = {}

    def register_reader(self, reader_name, reader_class):
        # for example, we can register the AlphaPept reader
        self.reader_dict[reader_name.lower()] = reader_class

    def get_reader(self, reader_name,
        frag_types=['b','y','b_modloss','y_modloss'],
        max_frag_charge=2,
        frag_tol=20, frag_ppm=True,
    )->PSMReader_w_FragBase:
        return self.reader_dict[reader_name.lower()](
            frag_types, max_frag_charge,
            frag_tol, frag_ppm,
        )

psm_w_frag_reader_provider = PSMwFragReaderProvider()