import torch
import pandas as pd
import numpy as np

import peptdeep.model.building_block as building_block
from peptdeep.model.model_interface import ModelInterface

ASCII_NUM=128

class Model_for_Generic_AASeq_Regression_LSTM(torch.nn.Module):
    """Generic LSTM regression model for AA sequence"""
    def __init__(self, 
        *,
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()
        self.dropout = torch.nn.Dropout(dropout)
        
        self.nn = torch.nn.Sequential(
            building_block.ascii_embedding(hidden_dim//4),
            building_block.SeqCNN(hidden_dim//4),
            self.dropout,
            building_block.SeqLSTM(
                hidden_dim, hidden_dim, 
                rnn_layer=nlayers
            ),
            building_block.SeqAttentionSum(hidden_dim),
            self.dropout,
            torch.nn.Linear(hidden_dim,64),
            torch.nn.GELU(),
            torch.nn.Linear(64, output_dim),
        )
    def forward(self, aa_x):
        return self.nn(aa_x).squeeze(-1)

class Model_for_Generic_AASeq_Regression_Transformer(torch.nn.Module):
    """Generic transformer regression model for AA sequence """
    def __init__(self,
        *,
        hidden_dim = 256,
        output_dim = 1,
        nlayers = 4,
        output_attentions=False,
        dropout = 0.1,
        **kwargs,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self.input_nn = building_block.ascii_embedding(hidden_dim)

        self.output_attentions = output_attentions
        
        self.hidden_nn = building_block.HFace_Transformer_with_PositionalEncoder(
            hidden_dim, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )

        self.output_nn = torch.nn.Sequential(
            building_block.SeqAttentionSum(hidden_dim),
            torch.nn.PReLU(),
            self.dropout,
            torch.nn.Linear(hidden_dim, output_dim),
        )

    @property
    def output_attentions(self)->bool:
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val:bool):
        self._output_attentions = val

    def forward(self, aa_x):
        aa_x = self.dropout(self.input_nn(aa_x))

        aa_x = self.hidden_nn(aa_x)
        if self.output_attentions:
            self.attentions = aa_x[1]
        else:
            self.attentions = None
        aa_x = self.dropout(aa_x[0])

        return self.output_nn(aa_x).squeeze(1)


class ModelInterface_for_Generic_AASeq_Regression(ModelInterface):
    """
    `ModelInterface` for Generic_AASeq_Regression models
    """
    def __init__(self, 
        model_class:torch.nn.Module=Model_for_Generic_AASeq_Regression_LSTM, 
        dropout=0.1,
        device:str='gpu',
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        **kwargs,
    ):
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            nlayers=nlayers,
            **kwargs
        )
        self.loss_func = torch.nn.L1Loss() # for regression

        self.target_column_to_predict = 'predicted_property'
        self.target_column_to_train = 'detected_property'

class Model_for_Generic_ModAASeq_Regression_LSTM(torch.nn.Module):
    """Generic LSTM regression model for modified sequence"""
    def __init__(self, 
        *,
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()
        self.dropout = torch.nn.Dropout(dropout)

        self.encoder_nn = building_block.Encoder_AA_Mod_CNN_LSTM_AttnSum(
            hidden_dim,
            n_lstm_layers=nlayers,
        )
        self.output_nn = torch.nn.Sequential(
            self.dropout,
            torch.nn.Linear(hidden_dim,64),
            torch.nn.GELU(),
            torch.nn.Linear(64, output_dim),
        )
    def forward(self, aa_x, mod_x):
        x = self.encoder_nn(aa_x, mod_x)
        return self.output_nn(x).squeeze(-1)

class Model_for_Generic_ModAASeq_Regression_Transformer(torch.nn.Module):
    """Generic transformer regression model for modified sequence"""
    def __init__(self,
        *,
        hidden_dim = 256,
        output_dim = 1,
        nlayers = 4,
        output_attentions=False,
        dropout = 0.1,
        **kwargs,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self.input_nn = building_block.AA_Mod_Embedding(hidden_dim)

        self._output_attentions = output_attentions
        
        self.hidden_nn = building_block.HFace_Transformer_with_PositionalEncoder(
            hidden_dim, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )

        self.output_nn = torch.nn.Sequential(
            building_block.SeqAttentionSum(hidden_dim),
            torch.nn.PReLU(),
            self.dropout,
            torch.nn.Linear(hidden_dim, output_dim),
        )

    @property
    def output_attentions(self)->bool:
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val:bool):
        self._output_attentions = val

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

class ModelInterface_for_Generic_ModAASeq_Regression(ModelInterface):
    """
    `ModelInterface` for all Generic_ModAASeq_Regression models
    """
    def __init__(self, 
        model_class:torch.nn.Module=Model_for_Generic_ModAASeq_Regression_LSTM,
        dropout=0.1,
        device:str='gpu',
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        **kwargs,
    ):
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            nlayers=nlayers,
            **kwargs
        )
        self.loss_func = torch.nn.L1Loss() # for regression

        self.target_column_to_predict = 'predicted_property'
        self.target_column_to_train = 'detected_property'

    def _get_features_from_batch_df(self, 
        batch_df: pd.DataFrame,
        **kwargs,
    ):
        return self._get_aa_mod_features(batch_df)

class Model_for_Generic_AASeq_BinaryClassification_LSTM(
    Model_for_Generic_AASeq_Regression_LSTM
):
    """Generic LSTM classification model for AA sequence"""
    def __init__(self, 
        *,
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__(
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            nlayers=nlayers,
            dropout=dropout,
        )

    def forward(self, aa_x):
        x = super().forward(aa_x)
        return torch.sigmoid(x)

class Model_for_Generic_AASeq_BinaryClassification_Transformer(
    Model_for_Generic_AASeq_Regression_Transformer
):
    """Generic transformer classification model for AA sequence"""
    def __init__(self,
        *,
        hidden_dim = 256,
        output_dim = 1,
        nlayers = 4,
        output_attentions=False,
        dropout = 0.1,
        **kwargs,
    ):
        """
        Model based on a transformer Architecture from 
        Huggingface's BertEncoder class.
        """
        super().__init__(
            nlayers=nlayers,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            output_attentions=output_attentions,
            dropout=dropout,
            **kwargs,
        )

    def forward(self, aa_x):
        x = super().forward(aa_x)
        return torch.sigmoid(x)

class ModelInterface_for_Generic_AASeq_BinaryClassification(ModelInterface):
    """
    `ModelInterface` for all Generic_AASeq_BinaryClassification models
    """
    def __init__(self, 
        model_class:torch.nn.Module=Model_for_Generic_AASeq_BinaryClassification_LSTM,
        dropout=0.1,
        device:str='gpu',
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        **kwargs,
    ):
        """
        Class to predict retention times from precursor dataframes.
        """
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            nlayers=nlayers,
            **kwargs
        )
        self.loss_func = torch.nn.BCELoss() # for binary classification
        self.target_column_to_predict = 'predicted_prob'
        self.target_column_to_train = 'detected_prob'

class Model_for_Generic_ModAASeq_BinaryClassification_LSTM(
    Model_for_Generic_ModAASeq_Regression_LSTM
):
    """Generic LSTM classification model for modified sequence"""
    def __init__(self, 
        *,
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__(
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            nlayers=nlayers,
            dropout=dropout,
            **kwargs,
        )

    def forward(self, aa_x, mod_x):
        x = super().forward(aa_x, mod_x)
        return torch.sigmoid(x)


class Model_for_Generic_ModAASeq_BinaryClassification_Transformer(
    Model_for_Generic_ModAASeq_Regression_Transformer
):
    """Generic transformer classification model for modified sequence"""
    def __init__(self,
        *,
        hidden_dim = 256,
        output_dim = 1,
        nlayers = 4,
        output_attentions=False,
        dropout = 0.1,
        **kwargs,
    ):
        super().__init__(
            nlayers=nlayers,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            output_attentions=output_attentions,
            dropout=dropout,
            **kwargs
        )

    @property
    def output_attentions(self)->bool:
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val:bool):
        self._output_attentions = val

    def forward(self, 
        aa_indices, 
        mod_x,
    ):
        x = super().forward(aa_indices, mod_x)
        return torch.sigmoid(x)


class ModelInterface_for_Generic_ModAASeq_BinaryClassification(ModelInterface):
    """
    `ModelInterface` for Generic_ModAASeq_BinaryClassification
    """
    def __init__(self, 
        model_class:torch.nn.Module=Model_for_Generic_ModAASeq_BinaryClassification_LSTM,
        dropout=0.1,
        device:str='gpu',
        hidden_dim=256,
        output_dim=1,
        nlayers=4,
        **kwargs,
    ):
        super().__init__(device=device)
        self.build(
            model_class,
            hidden_dim=hidden_dim,
            output_dim=output_dim,
            nlayers=nlayers,
            dropout=dropout,
            **kwargs
        )
        self.loss_func = torch.nn.BCELoss() # for classification

        self.target_column_to_predict = 'predicted_prob'
        self.target_column_to_train = 'detected_prob'

    def _get_features_from_batch_df(self, 
        batch_df: pd.DataFrame,
    ):
        return self._get_aa_mod_features(batch_df)
