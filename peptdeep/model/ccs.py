import torch
import pandas as pd
import numpy as np

from tqdm import tqdm

from alphabase.peptide.fragment import update_precursor_mz

from alphabase.peptide.mobility import (
    ccs_to_mobility_for_df,
    mobility_to_ccs_for_df
)

from peptdeep.model.featurize import (
    get_batch_aa_indices, 
    get_batch_mod_feature
)

from peptdeep.settings import model_const

import peptdeep.model.base as model_base


class Model_CCS_Bert(torch.nn.Module):
    """
    Transformer model for CCS prediction
    """
    def __init__(self,
        dropout = 0.1,
        nlayers = 4,
        hidden = 128,
        output_attentions=False,
        **kwargs,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self.input_nn = model_base.AATransformerEncoding(hidden-2)

        self._output_attentions = output_attentions
        
        self.hidden_nn = model_base.Hidden_HFace_Transformer(
            hidden, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )

        self.output_nn = torch.nn.Sequential(
            model_base.SeqAttentionSum(hidden),
            torch.nn.PReLU(),
            self.dropout,
            torch.nn.Linear(hidden, 1),
        )

    @property
    def output_attentions(self):
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val:bool):
        self._output_attentions = val
        self.hidden_nn.output_attentions = val

    def forward(self, 
        aa_indices, 
        mod_x,
        charges:torch.Tensor,
    ):
        x = self.dropout(self.input_nn(
            aa_indices, mod_x
        ))
        charges = charges.unsqueeze(1).repeat(1,x.size(1),2)
        x = torch.cat((x, charges),2)

        hidden_x = self.hidden_nn(x)
        if self.output_attentions:
            self.attentions = hidden_x[1]
        else:
            self.attentions = None
        x = self.dropout(hidden_x[0]+x*0.2)

        return self.output_nn(x).squeeze(1)

class Model_CCS_LSTM(torch.nn.Module):
    """LSTM model for CCS prediction"""
    def __init__(self,
        dropout=0.1
    ):
        super().__init__()
        
        self.dropout = torch.nn.Dropout(dropout)
        
        hidden = 256

        self.ccs_encoder = (
            model_base.Encoder_26AA_Mod_Charge_CNN_LSTM_AttnSum(
                hidden
            )
        )

        self.ccs_decoder = model_base.Decoder_Linear(
            hidden+1, 1
        )

    def forward(self, 
        aa_indices, 
        mod_x,
        charges,
    ):
        x = self.ccs_encoder(aa_indices, mod_x, charges)
        x = self.dropout(x)
        x = torch.cat((x, charges),1)
        return self.ccs_decoder(x).squeeze(1)

def ccs_to_mobility_pred_df(
    precursor_df:pd.DataFrame
)->pd.DataFrame:
    """ Add 'mobility_pred' into precursor_df inplace """
    precursor_df[
        'mobility_pred'
    ] = ccs_to_mobility_for_df(
        precursor_df, 'ccs_pred'
    )
    return precursor_df

def mobility_to_ccs_df_(
    precursor_df:pd.DataFrame
)->pd.DataFrame:
    """ Add 'ccs' into precursor_df inplace """
    precursor_df[
        'ccs'
    ] = mobility_to_ccs_for_df(
        precursor_df, 'mobility'
    )
    return precursor_df

class AlphaCCSModel(model_base.ModelInterface):
    """
    `ModelInterface` for `Model_CCS_LSTM` or `Model_CCS_Bert`
    """
    def __init__(self, 
        dropout=0.1,
        model_class:torch.nn.Module=Model_CCS_LSTM,
        device:str='gpu',
        **kwargs,
    ):
        super().__init__(device=device)
        self.model:Model_CCS_LSTM = None
        self.build(
            model_class,
            dropout=dropout, 
            **kwargs
        )
        self.charge_factor = 0.1

        self.target_column_to_predict = 'ccs_pred'
        self.target_column_to_train = 'ccs'

    def _get_features_from_batch_df(self, 
        batch_df: pd.DataFrame,
    ):
        aa_indices = self._get_26aa_indice_features(batch_df)

        mod_x = self._get_mod_features(batch_df)

        charges = self._as_tensor(
            batch_df['charge'].values
        ).unsqueeze(1)*self.charge_factor

        return aa_indices, mod_x, charges

    def ccs_to_mobility_pred(self,
        precursor_df:pd.DataFrame
    )->pd.DataFrame:
        return ccs_to_mobility_pred_df(precursor_df)
