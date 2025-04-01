import torch
import pandas as pd
import numpy as np
import warnings

from typing import List, Tuple, IO, Optional

from tqdm import tqdm

from alphabase.peptide.fragment import (
    init_fragment_by_precursor_dataframe,
    update_sliced_fragment_dataframe,
    get_sliced_fragment_dataframe,
    get_charged_frag_types,
    sort_charged_frag_types,
    parse_charged_frag_type,
    FRAGMENT_TYPES,
)

from peptdeep.utils import get_available_device

from peptdeep.model.featurize import (
    get_batch_aa_indices,
    parse_instrument_indices,
    get_batch_mod_feature,
)
from peptdeep.utils import logging
from peptdeep.settings import global_settings as settings, model_const

import peptdeep.model.model_interface as model_interface
import peptdeep.model.building_block as building_block


class ModelMS2Transformer(torch.nn.Module):
    """Transformer model for MS2 prediction

    Parameters
    ----------
    num_frag_types : int
        Total number of fragment types of a fragmentation position to predict

    num_modloss_types : int, optional
        Number of fragment types of a fragmentation position to predict, by default 0

    mask_modloss : bool, optional
        If True, the modloss layer will be disabled, by default True

    dropout : float, optional
        Dropout, by default 0.1

    nlayers : int, optional
        Number of transformer layer, by default 4

    hidden : int, optional
        Hidden layer size, by default 256
    """

    def __init__(
        self,
        num_frag_types: int,
        num_modloss_types: int = 0,
        mask_modloss: bool = True,
        dropout: float = 0.1,
        nlayers: int = 4,
        hidden: int = 256,
        **kwargs,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self._num_modloss_types = num_modloss_types
        self._num_non_modloss = num_frag_types - num_modloss_types
        self._mask_modloss = mask_modloss
        if num_modloss_types == 0:
            self._mask_modloss = True

        meta_dim = 8
        self.input_nn = building_block.Input_26AA_Mod_PositionalEncoding(
            hidden - meta_dim
        )

        self.meta_nn = building_block.Meta_Embedding(meta_dim)

        self.hidden_nn = building_block.Hidden_Transformer(
            hidden, nlayers=nlayers, dropout=dropout
        )

        self.output_nn = building_block.Decoder_Linear(
            hidden,
            self._num_non_modloss,
        )

        if num_modloss_types > 0:
            # for transfer learning of modloss frags
            self.modloss_nn = torch.nn.ModuleList(
                [
                    building_block.Hidden_Transformer(
                        hidden, nlayers=1, dropout=dropout
                    ),
                    building_block.Decoder_Linear(
                        hidden,
                        num_modloss_types,
                    ),
                ]
            )
        else:
            self.modloss_nn = None

    def forward(
        self,
        aa_indices,
        mod_x,
        charges: torch.Tensor,
        NCEs: torch.Tensor,
        instrument_indices,
    ):
        in_x = self.dropout(self.input_nn(aa_indices, mod_x))
        meta_x = (
            self.meta_nn(charges, NCEs, instrument_indices)
            .unsqueeze(1)
            .repeat(1, in_x.size(1), 1)
        )
        in_x = torch.cat((in_x, meta_x), 2)

        hidden_x = self.hidden_nn(in_x)
        hidden_x = self.dropout(hidden_x + in_x * 0.2)

        out_x = self.output_nn(hidden_x)

        if self._num_modloss_types > 0:
            if self._mask_modloss:
                out_x = torch.cat(
                    (
                        out_x,
                        torch.zeros(
                            *out_x.size()[:2],
                            self._num_modloss_types,
                            device=in_x.device,
                        ),
                    ),
                    2,
                )
            else:
                modloss_x = self.modloss_nn[0](in_x) + hidden_x
                modloss_x = self.modloss_nn[-1](modloss_x)
                out_x = torch.cat((out_x, modloss_x), 2)

        return out_x[:, 3:, :]


class ModelMS2Bert(torch.nn.Module):
    """Using HuggingFace's BertEncoder for MS2 prediction"""

    def __init__(
        self,
        charged_frag_types,
        dropout=0.1,
        nlayers=4,
        hidden=256,
        output_attentions=False,
        **kwargs,
    ):
        super().__init__()
        charged_frag_types = sort_charged_frag_types(charged_frag_types)
        self.dropout = torch.nn.Dropout(dropout)
        num_frag_types = len(charged_frag_types)

        # register charged fragment types
        self.register_buffer(
            "_supported_charged_frag_types",
            charged_frags_to_tensor(charged_frag_types),
        )
        self._get_modloss_frags()
        self._num_modloss_types = len(self._modloss_frag_types)
        self._num_non_modloss = num_frag_types - self._num_modloss_types

        meta_dim = 8
        self.input_nn = building_block.Input_26AA_Mod_PositionalEncoding(
            hidden - meta_dim
        )

        self.meta_nn = building_block.Meta_Embedding(meta_dim)

        self._output_attentions = output_attentions
        self.hidden_nn = building_block.Hidden_HFace_Transformer(
            hidden,
            nlayers=nlayers,
            dropout=dropout,
            output_attentions=output_attentions,
        )

        self.output_nn = building_block.Decoder_Linear(
            hidden,
            self._num_non_modloss,
        )

        if self._num_modloss_types > 0:
            # for transfer learning of modloss frags
            self.modloss_nn = torch.nn.ModuleList(
                [
                    building_block.Hidden_HFace_Transformer(
                        hidden,
                        nlayers=1,
                        dropout=dropout,
                        output_attentions=output_attentions,
                    ),
                    building_block.Decoder_Linear(
                        hidden,
                        self._num_modloss_types,
                    ),
                ]
            )
        else:
            self.modloss_nn = None

    def _get_modloss_frags(self):
        self._modloss_frag_types = []
        for i, frag in enumerate(self.supported_charged_frag_types):
            frag_type, _ = parse_charged_frag_type(frag)
            if FRAGMENT_TYPES[frag_type].modloss:
                self._modloss_frag_types.append(i)

    @property
    def output_attentions(self):
        return self._output_attentions

    @output_attentions.setter
    def output_attentions(self, val: bool):
        self._output_attentions = val
        self.hidden_nn.output_attentions = val
        self.modloss_nn[0].output_attentions = val

    @property
    def supported_charged_frag_types(self):
        return tensor_to_charged_frags(self._supported_charged_frag_types)

    def forward(
        self,
        aa_indices,
        mod_x,
        charges: torch.Tensor,
        NCEs: torch.Tensor,
        instrument_indices,
    ):
        in_x = self.dropout(self.input_nn(aa_indices, mod_x))
        meta_x = (
            self.meta_nn(charges, NCEs, instrument_indices)
            .unsqueeze(1)
            .repeat(1, in_x.size(1), 1)
        )
        in_x = torch.cat((in_x, meta_x), 2)

        hidden_x = self.hidden_nn(in_x)
        if self.output_attentions:
            self.attentions = hidden_x[1]
        else:
            self.attentions = None
        hidden_x = self.dropout(hidden_x[0] + in_x * 0.2)

        out_x = self.output_nn(hidden_x)

        self.modloss_attentions = None
        if self._num_modloss_types > 0:
            modloss_x = self.modloss_nn[0](in_x)
            if self.output_attentions:
                self.modloss_attentions = modloss_x[-1]
            modloss_x = modloss_x[0] + hidden_x
            modloss_x = self.modloss_nn[-1](modloss_x)
            out_x = torch.cat((out_x, modloss_x), 2)

        return out_x[:, 3:, :]


class ModelMS2pDeep(torch.nn.Module):
    """LSTM model for MS2 prediction similar to pDeep series"""

    def __init__(
        self,
        num_frag_types,
        num_modloss_types=0,
        mask_modloss=True,
        dropout=0.1,
        **kwargs,
    ):
        super().__init__()

        self.dropout = torch.nn.Dropout(dropout)

        self._num_modloss_types = num_modloss_types
        self._num_non_modloss = num_frag_types - num_modloss_types
        self._mask_modloss = mask_modloss
        if num_modloss_types == 0:
            self._mask_modloss = True

        BiRNN = True
        hidden = 512
        hidden_rnn_layer = 2

        self.input_nn = building_block.InputAALSTM_cat_Meta(hidden)

        self.hidden_nn = building_block.SeqLSTM(
            hidden, hidden, rnn_layer=hidden_rnn_layer, bidirectional=BiRNN
        )

        self.output_nn = building_block.OutputLSTM_cat_Meta(
            hidden,
            self._num_non_modloss,
        )

        if num_modloss_types:
            # for transfer learning of modloss frags
            self.modloss_nn = torch.nn.ModuleList(
                [
                    building_block.SeqLSTM(
                        hidden, hidden, rnn_layer=1, bidirectional=BiRNN
                    ),
                    building_block.SeqLSTM(
                        hidden, num_modloss_types, rnn_layer=1, bidirectional=False
                    ),
                ]
            )
        else:
            self.modloss_nn = None

    def forward(
        self,
        aa_indices,
        mod_x,
        charges: torch.Tensor,
        NCEs: torch.Tensor,
        instrument_indices,
    ):
        in_x = self.input_nn(aa_indices, mod_x, charges, NCEs, instrument_indices)
        in_x = self.dropout(in_x)

        hidden_x = self.hidden_nn(in_x)
        hidden_x = self.dropout(hidden_x)

        out_x = self.output_nn(hidden_x, charges, NCEs, instrument_indices)

        # modloss is mainly only for Phospho@S/T
        if self._num_modloss_types > 0:
            if self._mask_modloss:
                out_x = torch.cat(
                    (
                        out_x,
                        torch.zeros(
                            *out_x.size()[:2],
                            self._num_modloss_types,
                            device=in_x.device,
                        ),
                    ),
                    2,
                )
            else:
                modloss_x = self.modloss_nn[0](in_x) + hidden_x
                modloss_x = self.modloss_nn[-1](modloss_x)
                out_x = torch.cat((out_x, modloss_x), 2)

        return out_x[:, 3:, :]


class IntenAwareLoss(torch.nn.Module):
    """Loss weighted by intensity for MS2 models"""

    def __init__(self, base_weight=0.2):
        super().__init__()
        self.w = base_weight

    def forward(self, pred, target):
        return torch.mean((target + self.w) * torch.abs(target - pred))


mod_feature_size = len(model_const["mod_elements"])
max_instrument_num = model_const["max_instrument_num"]
frag_types = settings["model"]["frag_types"]
max_frag_charge = settings["model"]["max_frag_charge"]
num_ion_types = len(frag_types) * max_frag_charge


class pDeepModel(model_interface.ModelInterface):
    """
    `ModelInterface` for MS2 prediction models

    Parameters
    ----------
    charged_frag_types : List[str]
        Charged fragment types to predict
    dropout : float, optional
        Dropout rate, by default 0.1
    model_class : torch.nn.Module, optional
        Ms2 Model class, by default ModelMS2Bert
    device : str, optional
        Device to run the model, by default "gpu"
    override_from_weights : bool, optional default False
        Override the requested charged frag types from the model weights on loading. This allows to predict all fragment types supported by the weights even if the user doesn't know what fragments types are supported by the weights. Thereby, the model will always be in a safe to predict state.
    mask_modloss : bool, optional (deprecated)
        Mask the modloss fragments, this is deprecated and will be removed in the future. To mask the modloss fragments,
        the charged_frag_types should not include the modloss fragments.

    """

    def __init__(
        self,
        charged_frag_types=get_charged_frag_types(frag_types, max_frag_charge),
        dropout=0.1,
        model_class: torch.nn.Module = ModelMS2Bert,
        device: str = "gpu",
        mask_modloss: Optional[bool] = None,
        override_from_weights: bool = False,
        **kwargs,  # model params
    ):
        super().__init__(device=device)
        if mask_modloss is not None:
            warnings.warn(
                "mask_modloss is deprecated and will be removed in the future. To mask the modloss fragments, the charged_frag_types should not include the modloss fragments."
            )
        self.mask_modloss = mask_modloss
        self.override_from_weights = override_from_weights
        self.charged_frag_types = sort_charged_frag_types(charged_frag_types)

        self.charge_factor = 0.1
        self.NCE_factor = 0.01
        self.model: ModelMS2Bert = None
        self._model_kwargs = kwargs
        self.build(
            model_class,
            charged_frag_types=self.charged_frag_types,
            dropout=dropout,
            **kwargs,  # other model params
        )
        if mask_modloss:  # To Be removed in the future
            # remove modloss fragments from charged_frag_types
            self.charged_frag_types = [
                frag for frag in self.charged_frag_types if "modloss" not in frag
            ]

        self.loss_func = torch.nn.L1Loss()
        self.min_inten = 1e-4
        self._safe_to_predict = True
        self._safe_to_train = True

    def _prepare_train_data_df(
        self,
        precursor_df: pd.DataFrame,
        fragment_intensity_df: pd.DataFrame = None,
    ):
        self.frag_inten_df = fragment_intensity_df[self.charged_frag_types]
        # if np.all(precursor_df['nce'].values > 1):
        #     precursor_df['nce'] = precursor_df['nce']*self.NCE_factor

    def _check_predict_in_order(self, precursor_df: pd.DataFrame):
        pass

    def _prepare_predict_data_df(
        self,
        precursor_df: pd.DataFrame,
        reference_frag_df: pd.DataFrame = None,
    ):
        if reference_frag_df is None and precursor_df.nAA.is_monotonic_increasing:
            self._predict_in_order = True

            if "frag_start_idx" in precursor_df.columns:
                precursor_df.drop(
                    columns=["frag_start_idx", "frag_stop_idx"], inplace=True
                )
        else:
            self._predict_in_order = False

        self.predict_df = init_fragment_by_precursor_dataframe(
            precursor_df,
            self.charged_frag_types,
            reference_fragment_df=reference_frag_df,
            dtype=np.float32,
        )

        # if np.all(precursor_df['nce'].values > 1):
        #     precursor_df['nce'] = precursor_df['nce']*self.NCE_factor

    def _get_features_from_batch_df(
        self,
        batch_df: pd.DataFrame,
        **kwargs,
    ) -> Tuple[torch.Tensor]:
        aa_indices = self._get_26aa_indice_features(batch_df)

        mod_x = self._get_mod_features(batch_df)

        charges = (
            self._as_tensor(batch_df["charge"].values).unsqueeze(1) * self.charge_factor
        )

        nces = self._as_tensor(batch_df["nce"].values).unsqueeze(1) * self.NCE_factor

        instrument_indices = self._as_tensor(
            parse_instrument_indices(batch_df["instrument"]), dtype=torch.long
        )
        return aa_indices, mod_x, charges, nces, instrument_indices

    def _get_targets_from_batch_df(
        self, batch_df: pd.DataFrame, fragment_intensity_df: pd.DataFrame = None
    ) -> torch.Tensor:
        return self._as_tensor(
            get_sliced_fragment_dataframe(
                fragment_intensity_df,
                batch_df[["frag_start_idx", "frag_stop_idx"]].values,
            ).values
        ).view(-1, batch_df.nAA.values[0] - 1, len(self.charged_frag_types))

    def _set_batch_predict_data(
        self,
        batch_df: pd.DataFrame,
        predicts: np.ndarray,
        **kwargs,
    ):
        apex_intens = predicts.reshape((len(batch_df), -1)).max(axis=1)
        apex_intens[apex_intens <= 0] = 1
        predicts /= apex_intens.reshape((-1, 1, 1))
        predicts[predicts < self.min_inten] = 0.0
        # mask out predicted charged frag types that are not in the requested charged_frag_types
        columns_mask = np.isin(
            self.model.supported_charged_frag_types, self.charged_frag_types
        )
        predicts = predicts[:, :, columns_mask]

        if self._predict_in_order:
            self.predict_df.values[
                batch_df.frag_start_idx.values[0] : batch_df.frag_stop_idx.values[-1], :
            ] = predicts.reshape((-1, len(self.charged_frag_types)))
        else:
            update_sliced_fragment_dataframe(
                self.predict_df,
                self.predict_df.to_numpy(copy=True),
                predicts.reshape((-1, len(self.charged_frag_types))),
                batch_df[["frag_start_idx", "frag_stop_idx"]].values,
            )

    def _adapt_model_prediction_head(self):
        """
        Align the underlying model charged_frag_types with the interface charged_frag_types,
        this function is necessary for the model to be safe to train.
        Important: This function when called will reshape the prediction head of the model to match the requested charged_frag_types
        and randomly initialize is it so it might not be safe to use the model for prediction before training.
        """
        loaded_model_state_dict = self.model.state_dict()
        self.build(
            self.model.__class__,
            dropout=self.model.dropout.p,
            charged_frag_types=self.charged_frag_types,
            **self._model_kwargs,
        )
        current_model_dict = self.model.state_dict()
        # use the layers/modules that are shared from the original model instead of starting from scratch
        filtered_params = {}
        size_mismatches = []
        unexpected_keys = []
        for source_key, source_value in loaded_model_state_dict.items():
            if source_key in current_model_dict:
                if source_value.size() == current_model_dict[source_key].size():
                    filtered_params[source_key] = source_value
                else:
                    size_mismatches.append(source_key)
            else:
                unexpected_keys.append(source_key)
        self.model.load_state_dict(filtered_params, strict=False)

        if len(size_mismatches) > 0 or len(unexpected_keys) > 0:
            self._safe_to_predict = False

        self._safe_to_train = True

    def train_with_warmup(
        self,
        precursor_df: pd.DataFrame,
        fragment_intensity_df,
        *,
        batch_size=1024,
        epoch=10,
        warmup_epoch=5,
        lr=1e-5,
        verbose=False,
        verbose_each_epoch=False,
        **kwargs,
    ):
        if not self._safe_to_train:
            self._adapt_model_prediction_head()
        super().train_with_warmup(
            precursor_df,
            fragment_intensity_df=fragment_intensity_df,
            batch_size=batch_size,
            epoch=epoch,
            warmup_epoch=warmup_epoch,
            lr=lr,
            verbose=verbose,
            verbose_each_epoch=verbose_each_epoch,
            **kwargs,
        )
        self._safe_to_predict = True

    def test(
        self,
        precursor_df: pd.DataFrame,
        fragment_intensity_df: pd.DataFrame,
        default_instrument: str = "Lumos",
        default_nce: float = 30.0,
    ) -> pd.DataFrame:
        if "instrument" not in precursor_df.columns:
            precursor_df["instrument"] = default_instrument
        if "nce" not in precursor_df.columns:
            precursor_df["nce"] = default_nce
        columns = np.intersect1d(
            self.charged_frag_types,
            fragment_intensity_df.columns.values,
        )
        return calc_ms2_similarity(
            precursor_df,
            self.predict(precursor_df, reference_frag_df=fragment_intensity_df)[
                columns
            ],
            fragment_intensity_df=fragment_intensity_df[columns],
        )[-1]

    def train(
        self,
        precursor_df: pd.DataFrame,
        fragment_intensity_df,
        *,
        batch_size=1024,
        epoch=20,
        warmup_epoch=0,
        lr=1e-5,
        verbose=False,
        verbose_each_epoch=False,
        **kwargs,
    ):
        if not self._safe_to_train:
            self._adapt_model_prediction_head()
        super().train(
            precursor_df,
            fragment_intensity_df=fragment_intensity_df,
            batch_size=batch_size,
            epoch=epoch,
            warmup_epoch=warmup_epoch,
            lr=lr,
            verbose=verbose,
            verbose_each_epoch=verbose_each_epoch,
            **kwargs,
        )
        self._safe_to_predict = True

    def predict(
        self,
        precursor_df: pd.DataFrame,
        *,
        batch_size=1024,
        verbose=False,
        reference_frag_df=None,
        allow_unsafe_predictions=False,
        **kwargs,
    ) -> pd.DataFrame:
        """
        Predict MS2 fragment intensities

        Parameters
        ----------
        precursor_df : pd.DataFrame
            Precursor DataFrame
        batch_size : int, optional
            Batch size, by default 1024
        verbose : bool, optional
            Verbose, by default False
        reference_frag_df : pd.DataFrame, optional
            Reference fragment intensity DataFrame, by default None
        allow_unsafe_predictions : bool, optional
            Allow a newly a randomly initialized model to be used for prediction, by default False
        Returns
        -------
        pd.DataFrame
            Predicted fragment intensities
        """
        if not self._safe_to_predict and not allow_unsafe_predictions:
            raise ValueError(
                f"The model is not safe to use for prediction. This might mean that the requested charged_frag_types {self.charged_frag_types} are not a subset of the charged_frag_types used to train the loaded pretrained model {self.model.supported_charged_frag_types}. Please choose a subset of the supported charged_frag_types or retrain the model with the requested charged_frag_types."
            )
        return super().predict(
            precursor_df,
            batch_size=batch_size,
            verbose=verbose,
            reference_frag_df=reference_frag_df,
            **kwargs,
        )

    def predict_mp(self, **kwargs) -> pd.DataFrame:
        warnings.warn(
            "Please use pretrained_models.ModelManager::predict_all() "
            "for MS2 prediction with multiprocessing"
        )

    def bootstrap_nce_search(
        self,
        psm_df: pd.DataFrame,
        fragment_intensity_df: pd.DataFrame,
        nce_first=15,
        nce_last=45,
        nce_step=3,
        instrument="Lumos",
        charged_frag_types: List = None,
        metric="PCC>0.9",  # or 'median PCC'
        max_psm_subset=3000,
        n_bootstrap=3,
        callback=None,
    ):
        nce_list = []
        for i in range(n_bootstrap):
            nce, instrument = self.grid_nce_search(
                psm_df,
                fragment_intensity_df,
                nce_first,
                nce_last,
                nce_step,
                [instrument],
                charged_frag_types,
                metric,
                max_psm_subset,
                n_bootstrap,
                callback,
            )
            nce_list.append(nce)
        return np.median(nce_list), instrument

    def _load_model_from_stream(self, stream):
        to_be_loaded_state_dict = torch.load(stream, map_location=self.device)
        if "_supported_charged_frag_types" in to_be_loaded_state_dict:
            loaded_charged_frag_types = to_be_loaded_state_dict[
                "_supported_charged_frag_types"
            ]
            # build a model that has the same charged_frag_types as the loaded model
            self.build(
                self.model.__class__,
                dropout=self.model.dropout.p,
                charged_frag_types=tensor_to_charged_frags(loaded_charged_frag_types),
                **self._model_kwargs,
            )

            if self.override_from_weights:
                self.charged_frag_types = tensor_to_charged_frags(
                    loaded_charged_frag_types
                )

        (missing_keys, unexpect_keys) = self.model.load_state_dict(
            to_be_loaded_state_dict, strict=False
        )
        self._update_model_state()
        missing_keys = [
            key for key in missing_keys if "_supported_charged_frag_types" not in key
        ]
        if len(missing_keys) > 0:
            logging.warn(
                f"nn parameters {missing_keys} are MISSING while loading models in {self.__class__}"
            )
        if len(unexpect_keys) > 0:
            logging.warn(
                f"nn parameters {unexpect_keys} are UNEXPECTED while loading models in {self.__class__}"
            )

    def _update_model_state(self):
        """
        Update the model state "safe_to_predict" and "safe_to_train".
        Depending on the the interface (pDeepModel) charged frag types and he underlying model charged frag types,
        the model state will be updated.
        - safe_to_predict: True if the interface charged frag types are a subset of the underlying model charged frag types
        - safe_to_train: True if the interface charged frag types are the same as the underlying model charged frag types
        """
        self._safe_to_predict = set(self.charged_frag_types).issubset(
            set(self.model.supported_charged_frag_types)
        )
        self._safe_to_train = set(self.charged_frag_types) == set(
            self.model.supported_charged_frag_types
        )

    def grid_nce_search(
        self,
        psm_df: pd.DataFrame,
        fragment_intensity_df: pd.DataFrame,
        nce_first=15,
        nce_last=45,
        nce_step=3,
        search_instruments=["Lumos"],
        charged_frag_types: List = None,
        metric="PCC>0.9",  # or 'median PCC'
        max_psm_subset=1000000,
        callback=None,
    ):
        if len(psm_df) > max_psm_subset:
            psm_df = psm_df.sample(max_psm_subset).copy()
        best_pcc = -1
        best_nce = 0.0
        best_instrument = None
        if "median" in metric:
            metric_row = "50%"
        else:
            metric_row = ">0.90"
        search_instruments = set(
            [
                settings["model_mgr"]["instrument_group"][inst]
                for inst in search_instruments
            ]
        )
        for inst in search_instruments:
            for nce in np.arange(nce_first, nce_last + nce_step, nce_step):
                psm_df["nce"] = nce
                psm_df["instrument"] = inst
                predict_inten_df = self.predict(
                    psm_df, reference_frag_df=fragment_intensity_df
                )
                df, metrics = calc_ms2_similarity(
                    psm_df,
                    predict_inten_df,
                    fragment_intensity_df,
                    charged_frag_types=charged_frag_types,
                    metrics=["PCC"],
                )
                pcc = metrics.loc[metric_row, "PCC"]
                if pcc > best_pcc:
                    best_pcc = pcc
                    best_nce = nce
                    best_instrument = inst
        return best_nce, best_instrument


def normalize_fragment_intensities(
    psm_df: pd.DataFrame, frag_intensity_df: pd.DataFrame
):
    """Normalize the intensities to 0-1 values inplace

    Parameters
    ----------
    psm_df : pd.DataFrame
        PSM DataFrame

    frag_intensity_df : pd.DataFrame
        Fragment intensity DataFrame to be normalized.
        Intensities will be normalzied inplace.
    """
    frag_intensity_df_np = frag_intensity_df.to_numpy()
    for i, (frag_start_idx, frag_stop_idx) in enumerate(
        psm_df[["frag_start_idx", "frag_stop_idx"]].values
    ):
        intens = frag_intensity_df_np[frag_start_idx:frag_stop_idx]
        max_inten = np.max(intens)
        if max_inten > 0:
            intens /= max_inten
        frag_intensity_df_np[frag_start_idx:frag_stop_idx, :] = intens
    frag_intensity_df.loc[:] = frag_intensity_df_np


def pearson_correlation(x: torch.Tensor, y: torch.Tensor):
    """Compute pearson correlation between 2 batches of 1-D tensors

    Parameters
    ----------
    x : torch.Tensor
        Shape (Batch, n)

    y : torch.Tensor
        Shape (Batch, n)

    """
    return torch.cosine_similarity(
        x - x.mean(dim=1, keepdim=True), y - y.mean(dim=1, keepdim=True), dim=1
    )


# legacy
pearson = pearson_correlation


def spectral_angle(cos):
    cos[cos > 1] = 1
    return 1 - 2 * torch.arccos(cos) / np.pi


def _get_ranks(x: torch.Tensor, device) -> torch.Tensor:
    sorted_idx = x.argsort(dim=1)
    flat_idx = (
        sorted_idx + torch.arange(x.size(0), device=device).unsqueeze(1) * x.size(1)
    ).flatten()
    ranks = torch.zeros_like(flat_idx)
    ranks[flat_idx] = (
        torch.arange(x.size(1), device=device)
        .unsqueeze(0)
        .repeat(x.size(0), 1)
        .flatten()
    )
    ranks = ranks.reshape(x.size())
    ranks[x == 0] = 0
    return ranks


def spearman_correlation(x: torch.Tensor, y: torch.Tensor, device):
    """Compute spearman correlation between 2 batches of 1-D tensors

    Parameters
    ----------
    x : torch.Tensor
        Shape (Batch, n)

    y : torch.Tensor
        Shape (Batch, n)

    """
    x_rank = _get_ranks(x, device).to(torch.float32)
    y_rank = _get_ranks(y, device).to(torch.float32)

    n = x.size(1)
    upper = 6 * torch.sum((x_rank - y_rank).pow(2), dim=1)
    down = n * (n**2 - 1.0)
    return 1.0 - (upper / down)


# legacy
spearman = spearman_correlation


def add_cutoff_metric(metrics_describ, metrics_df, thres=0.9):
    vals = []
    for col in metrics_describ.columns.values:
        vals.append(
            metrics_df.loc[metrics_df[col] > thres, col].count() / len(metrics_df)
        )
    metrics_describ.loc[f">{thres:.2f}"] = vals
    return metrics_describ


def calc_ms2_similarity(
    psm_df: pd.DataFrame,
    predict_intensity_df: pd.DataFrame,
    fragment_intensity_df: pd.DataFrame,
    charged_frag_types: List = None,
    metrics=["PCC", "COS", "SA", "SPC"],
    GPU=True,
    batch_size=10240,
    verbose=False,
    spc_top_k=0,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    if GPU:
        device, _ = get_available_device()
    else:
        device = torch.device("cpu")

    if charged_frag_types is None or len(charged_frag_types) == 0:
        charged_frag_types = fragment_intensity_df.columns.values

    _grouped = psm_df.groupby("nAA")

    if verbose:
        batch_tqdm = tqdm(_grouped)
    else:
        batch_tqdm = _grouped

    for met in metrics:
        psm_df[met] = 0.0

    for nAA, df_group in batch_tqdm:
        for i in range(0, len(df_group), batch_size):
            batch_end = i + batch_size
            batch_df = df_group.iloc[i:batch_end, :]

            pred_intens = torch.tensor(
                get_sliced_fragment_dataframe(
                    predict_intensity_df,
                    batch_df[["frag_start_idx", "frag_stop_idx"]].values,
                    charged_frag_types,
                ).values,
                dtype=torch.float32,
                device=device,
            ).reshape(-1, (nAA - 1) * len(charged_frag_types))

            frag_intens = torch.tensor(
                get_sliced_fragment_dataframe(
                    fragment_intensity_df,
                    batch_df[["frag_start_idx", "frag_stop_idx"]].values,
                    charged_frag_types,
                ).values,
                dtype=torch.float32,
                device=device,
            ).reshape(-1, (nAA - 1) * len(charged_frag_types))

            if "PCC" in metrics:
                psm_df.loc[batch_df.index, "PCC"] = (
                    pearson_correlation(pred_intens, frag_intens).cpu().detach().numpy()
                )

            if "COS" in metrics or "SA" in metrics:
                cos = torch.cosine_similarity(pred_intens, frag_intens, dim=1)
                psm_df.loc[batch_df.index, "COS"] = cos.cpu().detach().numpy()

                if "SA" in metrics:
                    psm_df.loc[batch_df.index, "SA"] = (
                        spectral_angle(cos).cpu().detach().numpy()
                    )

            if "SPC" in metrics:
                if spc_top_k > 1 and spc_top_k < frag_intens.size(1):
                    sorted_idx = frag_intens.argsort(dim=1, descending=True)
                    flat_idx = (
                        sorted_idx[:, :spc_top_k]
                        + torch.arange(
                            frag_intens.size(0), dtype=torch.int, device=device
                        ).unsqueeze(1)
                        * frag_intens.size(1)
                    ).flatten()
                    pred_intens = pred_intens.flatten()[flat_idx].reshape(
                        sorted_idx.size(0), -1
                    )
                    frag_intens = frag_intens.flatten()[flat_idx].reshape(
                        sorted_idx.size(0), -1
                    )
                psm_df.loc[batch_df.index, "SPC"] = (
                    spearman_correlation(pred_intens, frag_intens, device)
                    .cpu()
                    .detach()
                    .numpy()
                )

    metrics_describ = psm_df[metrics].describe()
    add_cutoff_metric(metrics_describ, psm_df, thres=0.9)
    add_cutoff_metric(metrics_describ, psm_df, thres=0.75)

    torch.cuda.empty_cache()
    return psm_df, metrics_describ


def charged_frags_to_tensor(charged_frags: List[str]) -> torch.Tensor:
    """
    Convert a list of strings (charged fragment types, modloss fragment types) to a tensor

    Parameters
    ----------
    list : List[str]
        List of strings
    """
    seperator = ","
    string = seperator.join(charged_frags)
    return torch.tensor([ord(char) for char in string], dtype=torch.int32).unsqueeze(0)


def tensor_to_charged_frags(tensor: torch.Tensor) -> List[str]:
    """
    Convert a tensor to a list of strings (charged fragment types, modloss fragment types)

    Parameters
    ----------
    tensor : torch.Tensor
        Tensor of int32
    """
    string = "".join([chr(char) for char in tensor[0].tolist()])
    return string.split(",")
