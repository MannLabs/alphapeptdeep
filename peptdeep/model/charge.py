import pandas as pd
import numpy as np

from peptdeep.model.generic_property_prediction import (
    ModelInterface_for_Generic_AASeq_MultiLabelClassification,
    Model_for_Generic_AASeq_BinaryClassification_Transformer,
    ModelInterface_for_Generic_ModAASeq_MultiLabelClassification,
    Model_for_Generic_ModAASeq_BinaryClassification_Transformer,
)

MAX_SUPPORTED_CHARGE = 10
MIN_SUPPORTED_CHARGE = 1


class _ChargeModelInterface:
    def __init__(self, *args, **kwargs):
        raise TypeError("Interface class cannot be instantiated.")

    def predict_charges_as_prob(
        self,
        pep_df: pd.DataFrame,
        min_precursor_charge: int,
        max_precursor_charge: int,
    ):
        df = self.predict_mp(
            pep_df.copy(),
            batch_size=self.predict_batch_size,
        )
        df.rename(columns={"charge_probs": "charge_prob"}, inplace=True)
        df["charge"] = [
            self.charge_range[
                min_precursor_charge
                - self.min_supported_predict_charge : max_precursor_charge
                - self.min_supported_predict_charge
                + 1
            ]
        ] * len(df)
        df["charge_prob"] = df.charge_prob.apply(
            lambda x: x[
                min_precursor_charge
                - self.min_supported_predict_charge : max_precursor_charge
                - self.min_supported_predict_charge
                + 1
            ]
        )
        df = df.explode(["charge", "charge_prob"], ignore_index=True).dropna(
            subset=["charge"]
        )
        df["charge"] = df.charge.astype(np.int8)
        df["charge_prob"] = df.charge_prob.astype(np.float32)
        return df

    def predict_prob_for_charge(
        self,
        precursor_df: pd.DataFrame,
    ):
        if "charge" not in precursor_df.columns:
            raise KeyError("precursor_df must contain `charge` column")
        precursor_df = self.predict_mp(
            precursor_df,
            batch_size=self.predict_batch_size,
        )
        precursor_df["charge_prob"] = (
            precursor_df[["charge_probs", "charge"]]
            .apply(
                lambda x: x.iloc[0][x.iloc[1] - self.min_supported_predict_charge],
                axis=1,
            )
            .astype(np.float32)
        )
        precursor_df.drop(columns="charge_probs", inplace=True)
        return precursor_df

    def predict_and_clip_charges(
        self,
        pep_df: pd.DataFrame,
        min_precursor_charge: int,
        max_precursor_charge: int,
        charge_prob_cutoff: float,
    ):
        df = self.predict_mp(
            pep_df.copy(),
            batch_size=self.predict_batch_size,
        )
        df.rename(columns={"charge_probs": "charge_prob"}, inplace=True)
        df["charge"] = df.charge_prob.apply(
            lambda x: self.charge_range[x > charge_prob_cutoff]
        )
        df["charge_prob"] = df.charge_prob.apply(lambda x: x[x > charge_prob_cutoff])
        df = df.explode(["charge", "charge_prob"]).dropna(subset=["charge"])
        df["charge"] = df.charge.astype(np.int8)
        df = df.query(
            f"charge>={min_precursor_charge} and charge<={max_precursor_charge}"
        ).reset_index(drop=True)
        df["charge_prob"] = df.charge_prob.astype(np.float32)
        return df

    def create_charge_indicators(
        self, psm_df: pd.DataFrame, group_by_modseq: bool = True
    ):
        if group_by_modseq:
            grouped_psm_df = group_psm_df_by_modseq(
                psm_df,
                min_charge=self.min_supported_predict_charge,
                max_charge=self.max_supported_predict_charge,
            )
        else:
            grouped_psm_df = group_psm_df_by_sequence(
                psm_df,
                min_charge=self.min_supported_predict_charge,
                max_charge=self.max_supported_predict_charge,
            )
        return grouped_psm_df

    def test(
        self,
        precursor_df: pd.DataFrame,
        batch_size: int = 1024,
        min_precursor_charge: int = MIN_SUPPORTED_CHARGE,
        max_precursor_charge: int = MAX_SUPPORTED_CHARGE,
    ):
        return evaluate_charge_precision_recall(
            self.predict(precursor_df, batch_size=batch_size),
            charge_prob_cutoff=0.3,
            min_precursor_charge=min_precursor_charge,
            max_precursor_charge=max_precursor_charge,
        )


class ChargeModelForModAASeq(
    ModelInterface_for_Generic_ModAASeq_MultiLabelClassification, _ChargeModelInterface
):
    """
    ModelInterface for charge prediction for modified peptides

     Parameters
    ----------
    device : str, optional
        Device to use for training and prediction, by default "gpu"
    """

    def __init__(self, device: str = "gpu"):
        super().__init__(
            num_target_values=MAX_SUPPORTED_CHARGE - MIN_SUPPORTED_CHARGE + 1,
            model_class=Model_for_Generic_ModAASeq_BinaryClassification_Transformer,
            nlayers=4,
            hidden_dim=128,
            dropout=0.1,
            device=device,
        )

        self.target_column_to_predict = "charge_probs"
        self.target_column_to_train = "charge_indicators"
        self.min_supported_predict_charge = MIN_SUPPORTED_CHARGE
        self.max_supported_predict_charge = MAX_SUPPORTED_CHARGE
        self.charge_range = np.arange(
            MIN_SUPPORTED_CHARGE, MAX_SUPPORTED_CHARGE + 1, dtype=np.int8
        )
        self.predict_batch_size = 1024


class ChargeModelForAASeq(
    ModelInterface_for_Generic_AASeq_MultiLabelClassification, _ChargeModelInterface
):
    """
    ModelInterface for charge prediction for amino acid sequence

    Parameters
    ----------
    min_charge : int, optional
        Minimum charge to predict, by default 1
    max_charge : int, optional
        Maximum charge to predict, by default 6
    device : str, optional
        Device to use for training and prediction, by default "gpu"
    """

    def __init__(self, device: str = "gpu"):
        super().__init__(
            num_target_values=MAX_SUPPORTED_CHARGE - MIN_SUPPORTED_CHARGE + 1,
            model_class=Model_for_Generic_AASeq_BinaryClassification_Transformer,
            nlayers=4,
            hidden_dim=128,
            dropout=0.1,
            device=device,
        )

        self.target_column_to_predict = "charge_probs"
        self.target_column_to_train = "charge_indicators"
        self.min_supported_predict_charge = MIN_SUPPORTED_CHARGE
        self.max_supported_predict_charge = MAX_SUPPORTED_CHARGE
        self.charge_range = np.arange(
            MIN_SUPPORTED_CHARGE, MAX_SUPPORTED_CHARGE + 1, dtype=np.int8
        )
        self.predict_batch_size = 1024


def group_psm_df_by_sequence(
    psm_df: pd.DataFrame,
    min_charge: int,
    max_charge: int,
):
    return (
        psm_df.groupby("sequence")["charge"]
        .apply(
            lambda x: get_charge_indicators(
                set(x), min_charge=min_charge, max_charge=max_charge
            )
        )
        .reset_index(drop=False)
        .rename(columns={"charge": "charge_indicators"})
    )


def group_psm_df_by_modseq(
    psm_df: pd.DataFrame,
    min_charge: int,
    max_charge: int,
):
    return (
        psm_df.groupby(["sequence", "mods", "mod_sites"])["charge"]
        .apply(
            lambda x: get_charge_indicators(
                set(x), min_charge=min_charge, max_charge=max_charge
            )
        )
        .reset_index(drop=False)
        .rename(columns={"charge": "charge_indicators"})
    )


def get_charge_indicators(
    charge_list,
    min_charge: int,
    max_charge: int,
):
    charge_indicators = np.zeros(max_charge - min_charge + 1)
    for charge in charge_list:
        if charge <= max_charge and charge >= min_charge:
            charge_indicators[charge - min_charge] = 1.0
        else:
            raise ValueError(
                f"Charge {charge} is out of range [{min_charge}, {max_charge}]"
            )
    return charge_indicators


def evaluate_charge_precision_recall(
    df: pd.DataFrame,
    charge_prob_cutoff: float = 0.3,
    min_precursor_charge: int = 1,
    max_precursor_charge: int = 6,
):
    """
    Evaluate the precision and recall of charge prediction.

    Parameters
    ----------
    df : pd.DataFrame
        DataFrame containing the "charge_indicators" and "charge_probs" columns.
    charge_prob_cutoff : float, optional
        Cutoff for charge prediction, by default 0.3
    min_precursor_charge : int, optional
        Minimum charge to predict, by default 1
    max_precursor_charge : int, optional
        Maximum charge to predict, by default 10

    Returns
    -------
    pd.DataFrame
        DataFrame containing the precision and recall for each charge.
    """
    target_charges = df["charge_indicators"].apply(
        lambda x: x[
            min_precursor_charge - MIN_SUPPORTED_CHARGE : max_precursor_charge
            - MIN_SUPPORTED_CHARGE
            + 1
        ]
    )
    target_charges = np.array(target_charges.tolist())

    predicted_charges = df["charge_probs"].apply(
        lambda x: x[
            min_precursor_charge - MIN_SUPPORTED_CHARGE : max_precursor_charge
            - MIN_SUPPORTED_CHARGE
            + 1
        ]
    )
    predicted_charges = np.array(predicted_charges.tolist())
    # cutoff
    predicted_charges = predicted_charges >= charge_prob_cutoff

    matches = np.logical_and(target_charges, predicted_charges)
    precision = np.sum(matches) / (np.sum(predicted_charges) + 1e-10)
    recall = np.sum(matches) / (np.sum(target_charges) + 1e-10)

    results = {
        "precision": precision,
        "recall": recall,
    }

    return results
