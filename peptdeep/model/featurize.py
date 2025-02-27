import numpy as np
import pandas as pd
from typing import List, Union
import torch
from peptdeep.settings import (
    model_const,
    mod_feature_size,
    MOD_TO_FEATURE,
    mod_elements,
    mod_elem_to_idx,
    _parse_mod_formula,
    update_all_mod_features,
)
from rdkit import Chem

from transformers import AutoTokenizer, RobertaModel

from alphabase.smiles.peptide import PeptideSmilesEncoder
from typing import List, Optional, Union
from abc import ABC, abstractmethod
from torch.nn.utils.rnn import pad_sequence
def parse_mod_feature(
    nAA: int, mod_names: List[str], mod_sites: List[int]
) -> np.ndarray:
    """
    Get modification feature of a given peptide (len=nAA).
    Note that `site=0` is for peptide N-term modification,
    `site=1` is for peptide C-term modification, and
    `1<=site<=nAA` is for residue modifications on the peptide.

    Parameters
    ----------
    nAA : int
        the lenght of the peptide sequence

    mod_names : List[str]
        the modification names

    mod_sites : List[str]
        the modification sites corresponding
        to `mod_names` on the peptide

    Returns
    -------
    np.ndarray
        2-D feature array with shape `(nAA+2,mod_feature_size)`

    """
    mod_x = np.zeros((nAA + 2, mod_feature_size))
    if len(mod_names) > 0:
        for site, mod in zip(mod_sites, mod_names):
            mod_x[site] += MOD_TO_FEATURE[mod]
        # mod_x[mod_sites] = [MOD_TO_FEATURE[mod] for mod in mod_names]
    return mod_x


def get_batch_mod_feature(batch_df: pd.DataFrame) -> np.ndarray:
    """
    Parameters
    ----------
    batch_df : pd.DataFrame
        dataframe with 'sequence', 'mods', 'mod_sites' and 'nAA' columns.
        All sequence lengths must be the same, meaning that nAA values must be equal.

    Returns
    -------
    np.ndarray
        3-D tensor with shape (batch_size, nAA+2, mod_feature_size)
    """

    mod_features_list = batch_df.mods.str.split(";").apply(
        lambda mod_names: [MOD_TO_FEATURE[mod] for mod in mod_names if len(mod) > 0]
    )
    mod_sites_list = batch_df.mod_sites.str.split(";").apply(
        lambda mod_sites: [int(site) for site in mod_sites if len(site) > 0]
    )
    mod_x_batch = np.zeros(
        (len(batch_df), batch_df.nAA.values[0] + 2, mod_feature_size)
    )
    for i, (mod_feats, mod_sites) in enumerate(zip(mod_features_list, mod_sites_list)):
        if len(mod_sites) > 0:
            for site, feat in zip(mod_sites, mod_feats):
                # Process multiple mods on one site
                mod_x_batch[i, site, :] += feat
            # mod_x_batch[i,mod_sites,:] = mod_feats
    return mod_x_batch


def get_batch_aa_indices(seq_array: Union[List, np.ndarray]) -> np.ndarray:
    """
    Convert peptide sequences into AA ID array. ID=0 is reserved for masking,
    so ID of 'A' is 1, ID of 'B' is 2, ..., ID of 'Z' is 26 (maximum).
    Zeros are padded into the N- and C-term for each sequence.

    Parameters
    ----------
    seq_array : Union[List,np.ndarray]
        list or 1-D array of sequences with the same length

    Returns
    -------
    np.ndarray
        2-D `np.int32` array with the shape
        `(len(seq_array), len(seq_array[0])+2)`. Zeros is padded into the
        N- and C-term of each sequence, so the 1st-D is `len(seq_array[0])+2`.

    """
    x = np.array(seq_array).view(np.int32).reshape(len(seq_array), -1) - ord("A") + 1
    # padding zeros at the N- and C-term
    return np.pad(x, [(0, 0)] * (len(x.shape) - 1) + [(1, 1)])


def get_ascii_indices(seq_array: Union[List, np.ndarray]) -> np.ndarray:
    """
    Convert peptide sequences into ASCII code array.
    The values are from 0 to 127.
    Zeros are padded into the N- and C-term for each sequence.

    Parameters
    ----------
    seq_array : Union[List,np.ndarray]
        list or 1-D array of sequences.

    Returns
    -------
    np.ndarray
        2-D `np.int32` array with the shape
        `(len(seq_array), max seq length+2)`.
        For the the sequence whose length is shorter than max seq length,
        zeros are padded to the missing values.

    """

    x = np.array(seq_array).view(np.int32).reshape(len(seq_array), -1)
    return np.pad(x, [(0, 0)] * (len(x.shape) - 1) + [(1, 1)])


instrument_dict = dict(
    zip(
        [inst.upper() for inst in model_const["instruments"]],
        range(len(model_const["instruments"])),
    )
)
unknown_inst_index = model_const["max_instrument_num"] - 1


def parse_instrument_indices(instrument_list):
    instrument_list = [inst.upper() for inst in instrument_list]
    instrument_list = [inst for inst in instrument_list]
    return [
        instrument_dict[inst] if inst in instrument_dict else unknown_inst_index
        for inst in instrument_list
    ]

class FeatureExtractor(ABC):
    """
    Abstract class for feature extraction. All feature extractors should implement
    the get_batch_mod_feature, get_batch_aa_feature
    """
    
    @abstractmethod
    def get_batch_mod_feature(self, batch_df: pd.DataFrame)-> np.ndarray:
        pass

    @abstractmethod
    def get_batch_aa_feature(self, seq_array: Union[List, np.ndarray]) -> np.ndarray:
        pass    

class BertaFeatureExtractor(FeatureExtractor):
    def __init__(self):
        self.tokenizer = AutoTokenizer.from_pretrained("DeepChem/ChemBERTa-77M-MTR")
        self.model = RobertaModel.from_pretrained("DeepChem/ChemBERTa-77M-MTR")
        self.peptide_encoder = PeptideSmilesEncoder()
        self.padding_value = 0
    
    def _get_batch_features_from_smiles(self, smiles: List[List[str]]) -> np.ndarray:
        embeddings = []
        for smile in smiles:
            tokens = self.tokenizer(smile, return_tensors="pt", padding=True)
            with torch.no_grad():
                outputs = self.model(**tokens)
                embedding = outputs.last_hidden_state 
                embedding = embedding[:,0,:]
                embeddings.append(embedding)
                print(embedding.shape)

        # padding and batching the embeddings
        embeddings = pad_sequence(embeddings, batch_first=True, padding_value=self.padding_value)

        return embeddings
    

    def get_batch_aa_feature(self, seq_array: Union[List, np.ndarray]) -> np.ndarray:
        """
        Convert peptide sequences into BERTa features using only the amino acid information.
        """
        # get smiles for each sequence 
        smiles = [self.peptide_encoder.peptide_to_smiles_per_amino_acid(seq) for seq in seq_array]
        return self._get_batch_features_from_smiles(smiles)
    
    def get_batch_mod_feature(self, batch_df):
        """
        Get BERTa features of a given peptide sequence with modifications.
        """
        smiles = [self.peptide_encoder.peptide_to_smiles_per_amino_acid(seq, mods, mod_sites) for seq, mods, mod_sites in zip(batch_df.sequence, batch_df.mods, batch_df.mod_sites)]

        return self._get_batch_features_from_smiles(smiles)
