# AUTOGENERATED! DO NOT EDIT! File to edit: nbdev_nbs/model/generic_property_prediction.ipynb (unless otherwise specified).

__all__ = ['ASCII_NUM', 'Model_for_Generic_AASeq_Regression_LSTM', 'Model_for_Generic_AASeq_Regression_Transformer',
           'ModelInterface_for_Generic_AASeq_Regression', 'Model_for_Generic_ModAASeq_Regression_LSTM',
           'Model_for_Generic_ModAASeq_Regression_Transformer', 'ModelInterface_for_Generic_ModAASeq_Regression',
           'Model_for_Generic_AASeq_BinaryClassification_LSTM',
           'Model_for_Generic_AASeq_BinaryClassification_Transformer',
           'ModelInterface_for_Generic_AASeq_BinaryClassification',
           'Model_for_Generic_ModAASeq_BinaryClassification_LSTM',
           'Model_for_Generic_ModAASeq_BinaryClassification_Transformer',
           'ModelInterface_for_Generic_ModAASeq_BinaryClassification']

# Cell
import torch
import peptdeep.model.building_block as building_block
from peptdeep.model.model_interface import ModelInterface
import pandas as pd
import numpy as np

ASCII_NUM=128

# Cell

class Model_for_Generic_AASeq_Regression_LSTM(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
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
                rnn_layer=n_lstm_layers
            ),
            building_block.SeqAttentionSum(hidden_dim),
            self.dropout,
            torch.nn.Linear(hidden_dim,64),
            torch.nn.GELU(),
            torch.nn.Linear(64, 1),
        )
    def forward(self, aa_x):
        return self.nn(aa_x).squeeze(-1)

class Model_for_Generic_AASeq_Regression_Transformer(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim = 256,
        nlayers = 4,
        output_attentions=False,
        dropout = 0.1,
        **kwargs,
    ):
        """
        Model based on a transformer Architecture from
        Huggingface's BertEncoder class.
        """
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self.input_nn =  torch.nn.Sequential(
            building_block.ascii_embedding(hidden_dim),
        )

        self.output_attentions = output_attentions

        self.hidden_nn = building_block.HFace_Transformer_with_PositionalEncoder(
            hidden_dim, nlayers=nlayers, dropout=dropout,
            output_attentions=output_attentions
        )

        self.output_nn = torch.nn.Sequential(
            building_block.SeqAttentionSum(hidden_dim),
            torch.nn.PReLU(),
            self.dropout,
            torch.nn.Linear(hidden_dim, 1),
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
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=Model_for_Generic_AASeq_Regression_LSTM,
        device:str='gpu',
        **kwargs,
    ):
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
            **kwargs
        )
        self.loss_func = torch.nn.L1Loss() # for regression

        self.target_column_to_predict = 'predicted_property'
        self.target_column_to_train = 'detected_property'

# Cell
class Model_for_Generic_ModAASeq_Regression_LSTM(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()
        self.dropout = torch.nn.Dropout(dropout)

        self.encoder_nn = building_block.Encoder_AA_Mod_CNN_LSTM_AttnSum(
            hidden_dim,
            n_lstm_layers=n_lstm_layers,
        )
        self.output_nn = torch.nn.Sequential(
            self.dropout,
            torch.nn.Linear(hidden_dim,64),
            torch.nn.GELU(),
            torch.nn.Linear(64, 1),
        )
    def forward(self, aa_x, mod_x):
        x = self.encoder_nn(aa_x, mod_x)
        return self.output_nn(x).squeeze(-1)

class Model_for_Generic_ModAASeq_Regression_Transformer(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim = 256,
        nlayers = 4,
        output_attentions=False,
        dropout = 0.1,
        **kwargs,
    ):
        """
        Model based on a transformer Architecture from
        Huggingface's BertEncoder class.
        """
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
            torch.nn.Linear(hidden_dim, 1),
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
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=Model_for_Generic_ModAASeq_Regression_LSTM, #model defined above
        device:str='gpu',
        **kwargs,
    ):
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
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

# Cell
class Model_for_Generic_AASeq_BinaryClassification_LSTM(
    Model_for_Generic_AASeq_Regression_LSTM
):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__(
            hidden_dim=hidden_dim,
            n_lstm_layers=n_lstm_layers,
            dropout=dropout,
        )

    def forward(self, aa_x):
        x = super().forward(aa_x)
        return torch.sigmoid(x)

class Model_for_Generic_AASeq_BinaryClassification_Transformer(
    Model_for_Generic_AASeq_Regression_Transformer
):
    def __init__(self,
        *,
        hidden_dim = 256,
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
            output_attentions=output_attentions,
            dropout=dropout,
            **kwargs,
        )

    def forward(self, aa_x):
        x = super().forward(aa_x)
        return torch.sigmoid(x)

class ModelInterface_for_Generic_AASeq_BinaryClassification(ModelInterface):
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=Model_for_Generic_AASeq_BinaryClassification_LSTM, #model defined above
        device:str='gpu',
        **kwargs,
    ):
        """
        Class to predict retention times from precursor dataframes.
        """
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
            **kwargs
        )
        self.loss_func = torch.nn.BCELoss() # for binary classification
        self.target_column_to_predict = 'predicted_prob'
        self.target_column_to_train = 'detected_prob'

# Cell
class Model_for_Generic_ModAASeq_BinaryClassification_LSTM(
    Model_for_Generic_ModAASeq_Regression_LSTM
):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__(
            hidden_dim=hidden_dim,
            n_lstm_layers=n_lstm_layers,
            dropout=dropout,
            **kwargs,
        )

    def forward(self, aa_x, mod_x):
        x = super().forward(aa_x, mod_x)
        return torch.sigmoid(x)


class Model_for_Generic_ModAASeq_BinaryClassification_Transformer(
    Model_for_Generic_ModAASeq_Regression_Transformer
):
    def __init__(self,
        *,
        hidden_dim = 256,
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
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=Model_for_Generic_ModAASeq_BinaryClassification_LSTM, #model defined above
        device:str='gpu',
        **kwargs,
    ):
        super().__init__(device=device)
        self.build(
            model_class,
            dropout=dropout,
            **kwargs
        )
        self.loss_func = torch.nn.BCELoss() # for regression

        self.target_column_to_predict = 'predicted_prob'
        self.target_column_to_train = 'detected_prob'

    def _get_features_from_batch_df(self,
        batch_df: pd.DataFrame,
    ):
        return self._get_aa_mod_features(batch_df)