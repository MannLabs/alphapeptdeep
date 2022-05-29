# AUTOGENERATED! DO NOT EDIT! File to edit: nbdev_nbs/model/model_shop.ipynb (unless otherwise specified).

__all__ = ['ASCII_NUM', 'ScalarRegression_LSTM_Model_for_AASeq', 'ScalarRegression_Transformer_Model_for_AASeq',
           'ScalarRegression_ModelInterface_for_AASeq', 'BinaryClassification_LSTM_Model_for_AASeq',
           'BinaryClassification_Transformer_Model_for_AASeq', 'BinaryClassification_ModelInterface_for_AASeq',
           'ScalarRegression_LSTM_Model_for_ModAASeq', 'ScalarRegression_Transformer_Model_for_ModAASeq',
           'ScalarRegression_ModelInterface_for_ModAASeq', 'BinaryClassification_LSTM_Model_for_ModAASeq',
           'BinaryClassification_Transformer_Model_for_ModAASeq', 'BinaryClassification_ModelInterface_for_ModAASeq']

# Cell
import torch
import peptdeep.model.building_block as building_block
from peptdeep.model.model_interface import ModelInterface
from peptdeep.model.featurize import (
    get_ascii_indices, get_batch_mod_feature
)
import pandas as pd
import numpy as np

ASCII_NUM=128

# Cell

class ScalarRegression_LSTM_Model_for_AASeq(torch.nn.Module):
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

class ScalarRegression_Transformer_Model_for_AASeq(torch.nn.Module):
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
    def output_attentions(self):
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val:bool):
        self._output_attentions = val
        self.hidden_nn.output_attentions = val

    def forward(self, aa_x):
        aa_x = self.dropout(self.input_nn(aa_x))

        aa_x = self.hidden_nn(aa_x)
        if self.output_attentions:
            self.attentions = aa_x[1]
        else:
            self.attentions = None
        aa_x = self.dropout(aa_x[0])

        return self.output_nn(aa_x).squeeze(1)

class ScalarRegression_ModelInterface_for_AASeq(ModelInterface):
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=ScalarRegression_LSTM_Model_for_AASeq, #model defined above
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

    def _prepare_predict_data_df(self,
        precursor_df:pd.DataFrame,
    ):
        self._predict_column_in_df = 'target_value_pred'
        precursor_df[self._predict_column_in_df] = 0.
        self.predict_df = precursor_df

    def _get_features_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs,
    ):
        aa_indices = self._as_tensor(
            get_ascii_indices(
                batch_df['sequence'].values.astype('U')
            ),
            dtype=torch.long
        )

        return aa_indices

    def _get_targets_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs
    ) -> torch.Tensor:
        return self._as_tensor(
            batch_df['target_value'].values,
            dtype=torch.float32
        )

# Cell
class BinaryClassification_LSTM_Model_for_AASeq(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()

        self.nn = ScalarRegression_LSTM_Model_for_AASeq(
            hidden_dim=hidden_dim,
            input_dim=ASCII_NUM,
            n_lstm_layers=n_lstm_layers,
            dropout=dropout,
        )
    def forward(self, aa_x):
        return torch.sigmoid(self.nn(aa_x))

class BinaryClassification_Transformer_Model_for_AASeq(torch.nn.Module):
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

        self.nn =  ScalarRegression_Transformer_Model_for_AASeq(
            nlayers=nlayers,
            input_dim=ASCII_NUM,
            hidden_dim=hidden_dim,
            output_attentions=output_attentions,
            dropout=dropout,
            **kwargs,
        )

    @property
    def output_attentions(self):
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val:bool):
        self._output_attentions = val
        self.hidden_nn.output_attentions = val

    def forward(self, aa_x):
        return torch.sigmoid(self.nn(aa_x))

class BinaryClassification_ModelInterface_for_AASeq(ModelInterface):
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=BinaryClassification_LSTM_Model_for_AASeq, #model defined above
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

    def _prepare_predict_data_df(self,
        precursor_df:pd.DataFrame,
    ):
        self._predict_column_in_df = 'target_prob_pred'
        precursor_df[self._predict_column_in_df] = 0.
        self.predict_df = precursor_df

    def _get_features_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs,
    ):
        aa_indices = self._as_tensor(
            get_ascii_indices(
                batch_df['sequence'].values.astype('U')
            ), dtype=torch.long
        )

        return aa_indices

    def _get_targets_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs
    ) -> torch.Tensor:
        return self._as_tensor(
            batch_df['target_prob'].values,
            dtype=torch.float32
        )

# Cell
class ScalarRegression_LSTM_Model_for_ModAASeq(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()
        self.dropout = torch.nn.Dropout(dropout)

        self.encoder_nn = building_block.Encoder_AsciiAA_Mod_CNN_LSTM_AttnSum(
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

class ScalarRegression_Transformer_Model_for_ModAASeq(torch.nn.Module):
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

class ScalarRegression_ModelInterface_for_ModAASeq(ModelInterface):
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=ScalarRegression_LSTM_Model_for_ModAASeq, #model defined above
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

    def _prepare_predict_data_df(self,
        precursor_df:pd.DataFrame,
    ):
        self._predict_column_in_df = 'target_value_pred'
        precursor_df[self._predict_column_in_df] = 0.
        self.predict_df = precursor_df

    def _get_features_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs,
    ):
        aa_indices = self._as_tensor(
            get_ascii_indices(
                batch_df['sequence'].values.astype('U')
            ),
            dtype=torch.long
        )

        mod_x = self._as_tensor(
            get_batch_mod_feature(
                batch_df
            )
        )

        return aa_indices, mod_x

    def _get_targets_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs
    ) -> torch.Tensor:
        return self._as_tensor(
            batch_df['target_value'].values,
            dtype=torch.float32
        )

# Cell
class BinaryClassification_LSTM_Model_for_ModAASeq(torch.nn.Module):
    def __init__(self,
        *,
        hidden_dim=256,
        n_lstm_layers=4,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()
        self.nn = ScalarRegression_LSTM_Model_for_ModAASeq(
            hidden_dim=hidden_dim,
            n_lstm_layers=n_lstm_layers,
            dropout=dropout
        )

    def forward(self, aa_x, mod_x):
        return torch.sigmoid(self.nn(aa_x, mod_x))

class BinaryClassification_Transformer_Model_for_ModAASeq(torch.nn.Module):
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
        self.nn = ScalarRegression_Transformer_Model_for_ModAASeq(
            nlayers=nlayers,
            hidden_dim=hidden_dim,
            output_attentions=output_attentions,
            dropout=dropout,
            **kwargs
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
        return torch.sigmoid(self.nn(aa_indices, mod_x))

class BinaryClassification_ModelInterface_for_ModAASeq(ModelInterface):
    def __init__(self,
        dropout=0.1,
        model_class:torch.nn.Module=BinaryClassification_LSTM_Model_for_ModAASeq, #model defined above
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

    def _prepare_predict_data_df(self,
        precursor_df:pd.DataFrame,
    ):
        self._predict_column_in_df = 'target_prob_pred'
        precursor_df[self._predict_column_in_df] = 0.
        self.predict_df = precursor_df

    def _get_features_from_batch_df(self,
        batch_df: pd.DataFrame,
    ):
        aa_indices = self._as_tensor(
            get_ascii_indices(
                batch_df['sequence'].values.astype('U')
            ),
            dtype=torch.long
        )

        mod_x = self._as_tensor(
            get_batch_mod_feature(
                batch_df
            )
        )

        return aa_indices, mod_x

    def _get_targets_from_batch_df(self,
        batch_df: pd.DataFrame,
        **kwargs
    ) -> torch.Tensor:
        return self._as_tensor(
            batch_df['target_prob'].values,
            dtype=torch.float32
        )