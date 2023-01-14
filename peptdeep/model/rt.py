import torch
import pandas as pd
import numpy as np

from peptdeep.model.featurize import (
    get_batch_aa_indices, 
    get_batch_mod_feature
)

from peptdeep.settings import model_const

import peptdeep.model.model_interface as model_interface
import peptdeep.model.building_block as building_block

mod_feature_size = len(model_const['mod_elements'])

IRT_PEPTIDE_DF = pd.DataFrame(
    [['LGGNEQVTR', 'RT-pep a', -24.92, '', ''],
    ['GAGSSEPVTGLDAK', 'RT-pep b', 0.00, '', ''],
    ['VEATFGVDESNAK', 'RT-pep c', 12.39, '', ''],
    ['YILAGVENSK', 'RT-pep d', 19.79, '', ''],
    ['TPVISGGPYEYR', 'RT-pep e', 28.71, '', ''],
    ['TPVITGAPYEYR', 'RT-pep f', 33.38, '', ''],
    ['DGLDAASYYAPVR', 'RT-pep g', 42.26, '', ''],
    ['ADVTPADFSEWSK', 'RT-pep h', 54.62, '', ''],
    ['GTFIIDPGGVIR', 'RT-pep i', 70.52, '', ''],
    ['GTFIIDPAAVIR', 'RT-pep k', 87.23, '', ''],
    ['LFLQFGAQGSPFLK', 'RT-pep l', 100.00, '', '']],
    columns=['sequence','pep_name','irt', 'mods', 'mod_sites']
)
IRT_PEPTIDE_DF['nAA'] = IRT_PEPTIDE_DF.sequence.str.len()


#legacy
irt_pep = IRT_PEPTIDE_DF


class Model_RT_Bert(torch.nn.Module):
    """Transformer model for RT prediction"""
    def __init__(self,
        dropout = 0.1,
        nlayers = 4,
        hidden = 128,
        output_attentions=False,
        **kwargs,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self.input_nn = building_block.AATransformerEncoding(hidden)

        self._output_attentions = output_attentions
        
        self.hidden_nn = building_block.Hidden_HFace_Transformer(
            hidden, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )

        self.output_nn = torch.nn.Sequential(
            building_block.SeqAttentionSum(hidden),
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
    ):
        x = self.dropout(self.input_nn(
            aa_indices, mod_x
        ))

        hidden_x = self.hidden_nn(x)
        if self.output_attentions:
            self.attentions = hidden_x[1]
        else:
            self.attentions = None
        x = self.dropout(hidden_x[0]+x*0.2)

        return self.output_nn(x).squeeze(1)


class Model_RT_LSTM_CNN(torch.nn.Module):
    """CNN+LSTM model for RT prediction"""
    def __init__(self, 
        dropout=0.2,
    ):
        super().__init__()
        
        self.dropout = torch.nn.Dropout(dropout)
        
        hidden = 256
        self.rt_encoder = building_block.Encoder_26AA_Mod_CNN_LSTM_AttnSum(
            hidden
        )

        self.rt_decoder = building_block.Decoder_Linear(
            hidden,
            1
        )

    def forward(self, 
        aa_indices, 
        mod_x,
    ):
        x = self.rt_encoder(aa_indices, mod_x)
        x = self.dropout(x)

        return self.rt_decoder(x).squeeze(1)
#legacy
Model_RT_LSTM = Model_RT_LSTM_CNN

class AlphaRTModel(model_interface.ModelInterface):
    """
    `ModelInterface` for RT models
    """
    def __init__(self, 
        dropout=0.1,
        model_class:torch.nn.Module=Model_RT_LSTM_CNN, #model defined above
        device:str='gpu',
        **kwargs,
    ):
        super().__init__(device=device)
        self.model:Model_RT_LSTM_CNN = None
        self.build(
            model_class,
            dropout=dropout,
            **kwargs
        )
        self.target_column_to_predict = 'rt_pred'
        self.target_column_to_train = 'rt_norm'

    def _get_features_from_batch_df(self, 
        batch_df: pd.DataFrame,
    ):
        return (
            self._get_26aa_indice_features(batch_df),
            self._get_mod_features(batch_df)
        )

    def add_irt_column_to_precursor_df(self,
        precursor_df: pd.DataFrame
    ):
        self.predict(IRT_PEPTIDE_DF)
        # simple linear regression
        rt_pred_mean = IRT_PEPTIDE_DF.rt_pred.mean()
        irt_mean = IRT_PEPTIDE_DF.irt.mean()
        x = IRT_PEPTIDE_DF.rt_pred.values - rt_pred_mean
        y = IRT_PEPTIDE_DF.irt.values - irt_mean
        slope = np.sum(x*y)/np.sum(x*x)
        intercept = irt_mean - slope*rt_pred_mean
        # end linear regression
        precursor_df['irt_pred'] = precursor_df.rt_pred*slope + intercept
        return precursor_df


