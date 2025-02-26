import pytest
from peptdeep.model.ms2 import pDeepModel
import os
import tempfile
from typing import List, Optional
from peptdeep.model.rt import IRT_PEPTIDE_DF
import numpy as np


def get_legacy_model(charged_frag_types: Optional[List[str]] = None):
    if charged_frag_types is None:
        return pDeepModel()
    else:
        return pDeepModel(charged_frag_types=charged_frag_types)


def transform_weights_to_new_format():
    # TODO This is a temporary solution to transform the old weights to the new format until the new weights are uploaded
    # load the legacy model with the correct charged fragment types and resave them to the new format
    model = get_legacy_model(
        charged_frag_types=[
            "b_z1",
            "b_z2",
            "y_z1",
            "y_z2",
            "b_modloss_z1",
            "b_modloss_z2",
            "y_modloss_z1",
            "y_modloss_z2",
        ]
    )
    temp_dir = os.path.join(tempfile.gettempdir(), "peptDeep_models")
    os.makedirs(temp_dir, exist_ok=True)
    weights_dist = os.path.join(temp_dir, "new_weights.pth")
    model.save(weights_dist)

    return weights_dist


def get_prediction_dataset():
    df = IRT_PEPTIDE_DF.copy()
    df["charge"] = 2
    df["mods"] = ""
    df["mod_sites"] = ""
    # sort by nAA
    df = df.sort_values("nAA")
    idxes = np.zeros(len(df) + 1, dtype=np.int64)
    idxes[1:] = np.cumsum(df.nAA.values - 1)
    df["frag_start_idx"] = idxes[:-1]
    df["frag_stop_idx"] = idxes[1:]
    df["nce"] = 30
    df["instrument"] = "Lumos"
    # sort by
    return df


def test_legacy_weights_with_correct_frag_types():
    # Given a user requests exactly the same charged frag types used when training the model
    charged_frag_types = [
        "b_z1",
        "b_z2",
        "y_z1",
        "y_z2",
        "b_modloss_z1",
        "b_modloss_z2",
        "y_modloss_z1",
        "y_modloss_z2",
    ]

    # When the user loads the model from legacy weights
    model = get_legacy_model(charged_frag_types=charged_frag_types)

    # Then the model should be safe to predict and train
    assert model._safe_to_predict, (
        "Model was not safe to predict when loading legacy weights with correct charged frag types"
    )
    assert model._safe_to_train, (
        "Model was not safe to train when loading legacy weights with correct charged frag types"
    )


def test_legacy_weights_complete_prediction():
    # Given a user requests exactly the same charged frag types used when training the model
    charged_frag_types = [
        "b_z1",
        "b_z2",
        "y_z1",
        "y_z2",
        "b_modloss_z1",
        "b_modloss_z2",
        "y_modloss_z1",
        "y_modloss_z2",
    ]

    # When the user loads the model from legacy weights and uses it for prediction
    model = get_legacy_model(charged_frag_types=charged_frag_types)
    df = get_prediction_dataset()
    pred_df = model.predict(df)

    # Then the prediction should have all requested charged frag types
    assert set(pred_df.columns) == set(charged_frag_types), (
        "Prediction did not have all requested charged frag types"
    )
    # Non nan values should be present
    assert not pred_df.isna().all().all(), "All values in the prediction were nan"


def test_new_weights_complete_prediction():
    # Given a user requests exactly the same charged frag types used when training the model
    requested_charged_frag_types = [
        "b_z1",
        "b_z2",
        "y_z1",
        "y_z2",
        "b_modloss_z1",
        "b_modloss_z2",
        "y_modloss_z1",
        "y_modloss_z2",
    ]

    # When the user loads the model from new weights and uses it for prediction
    model = pDeepModel(charged_frag_types=requested_charged_frag_types)
    model.load(transform_weights_to_new_format())

    df = get_prediction_dataset()

    pred_df = model.predict(df)

    # Then the prediction should have all requested charged frag types
    assert set(pred_df.columns) == set(requested_charged_frag_types), (
        "Prediction did not have all requested charged frag types"
    )
    # Non nan values should be present
    assert not pred_df.isna().all().all(), "All values in the prediction were nan"


def test_new_state_subset_prediction():
    # Given a user requests a subset of the charged frag types used when training the model
    requested_charged_frag_types = ["b_z1", "b_z2", "y_z1", "y_z2"]

    # When the user loads the model from new weights
    model = pDeepModel(charged_frag_types=requested_charged_frag_types)
    model.load(transform_weights_to_new_format())

    # The model should be safe to predict but not to train
    assert model._safe_to_predict, (
        "Model was not safe to predict when loading new weights with subset of charged frag types"
    )
    # since theres a discrepancy between the requested and the trained charged frag types
    assert not model._safe_to_train, (
        "Model was safe to train when loading new weights with subset of charged frag types"
    )


def test_new_state_unsupported_frag_types():
    # Given a user requests a charged frag types that are not supported by the loaded model
    requested_charged_frag_types = [
        "b_z1",
        "b_z2",
        "y_z1",
        "y_z2",
        "b_modloss_z1",
        "b_modloss_z2",
        "y_modloss_z1",
        "y_modloss_z2",
        "b_modloss_z3",
    ]

    # When the user loads the model from new weights
    model = pDeepModel(charged_frag_types=requested_charged_frag_types)
    model.load(transform_weights_to_new_format())

    # The model should not be safe to predict or train
    assert not model._safe_to_predict, (
        "Model was safe to predict when loading new weights with unsupported charged frag types"
    )
    assert not model._safe_to_train, (
        "Model was safe to train when loading new weights with unsupported charged frag types"
    )


def test_prediction_unsupported_frag_types():
    # Given a user requests a charged frag types that are not supported by the loaded model
    requested_charged_frag_types = ["b_z1", "b_z2", "x_z1", "x_z2"]

    # When the user loads the model from new weights and uses it for prediction
    model = pDeepModel(requested_charged_frag_types)
    model.load(transform_weights_to_new_format())

    df = get_prediction_dataset()
    # Then the model should raise an error
    with pytest.raises(ValueError):
        model.predict(df)


def test_override_requested_frag_types():
    # Given a user requests any charged frag types
    requested_charged_frag_types = ["b_z1", "b_z2", "x_z1", "x_z2"]

    # When the user loads the model from new weights and uses it for prediction while using the supported charged frag types to override the requested ones
    model = pDeepModel(requested_charged_frag_types, override_from_weights=True)
    model.load(transform_weights_to_new_format())

    # Then the requested charged frag types should be overridden by the supported ones and the model should be safe to predict
    assert set(model.charged_frag_types) == set(
        [
            "b_z1",
            "b_z2",
            "y_z1",
            "y_z2",
            "b_modloss_z1",
            "b_modloss_z2",
            "y_modloss_z1",
            "y_modloss_z2",
        ]
    ), "Overridden charged frag types were not as expected"
    assert model._safe_to_predict, (
        "Model was not safe to predict when overriding requested charged frag types"
    )


def test_override_requested_frag_types_prediction():
    # Given a user requests any charged frag types
    requested_charged_frag_types = ["b_z1", "b_z2", "x_z1", "x_z2"]

    # When the user loads the model from new weights and uses it for prediction while using the supported charged frag types to override the requested ones
    model = pDeepModel(requested_charged_frag_types, override_from_weights=True)
    model.load(transform_weights_to_new_format())
    df = get_prediction_dataset()
    pred_df = model.predict(df)

    # Then the prediction should have all charged frag types supported by the model
    assert set(pred_df.columns) == set(
        [
            "b_z1",
            "b_z2",
            "y_z1",
            "y_z2",
            "b_modloss_z1",
            "b_modloss_z2",
            "y_modloss_z1",
            "y_modloss_z2",
        ]
    ), "Prediction did not have all supported charged frag types"
    # Non nan values should be present
    assert not pred_df.isna().all().all(), "All values in the prediction were nan"
