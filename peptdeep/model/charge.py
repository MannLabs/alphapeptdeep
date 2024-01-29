import pandas as pd
import numpy as np


from peptdeep.model.generic_property_prediction import (
    ModelInterface_for_Generic_AASeq_MultiLabelClassification,
    Model_for_Generic_AASeq_BinaryClassification_Transformer,
    ModelInterface_for_Generic_ModAASeq_MultiLabelClassification,
    Model_for_Generic_ModAASeq_BinaryClassification_Transformer,
)
    
class ChargeModelForModAASeq(
    ModelInterface_for_Generic_ModAASeq_MultiLabelClassification
):
    """
    ModelInterface for charge prediction for modified peptides
    """
    def __init__(self, min_charge:int=1, max_charge:int=6):
        super().__init__(
            num_target_values=max_charge-min_charge+1,
            model_class=Model_for_Generic_ModAASeq_BinaryClassification_Transformer,
            nlayers=4, hidden_dim=128, dropout=0.1
        )

        self.target_column_to_predict = "charge_probs"
        self.target_column_to_train = "charge_indicators"
        self.min_charge = min_charge
        self.max_charge = max_charge
        self.charge_range = np.arange(
            min_charge, max_charge+1, dtype=np.int8
        )
        
    def predict_charges_for_pep_df(self, 
        pep_df:pd.DataFrame, 
        charge_prob=0.3,
        drop_probs_column=True
    ):
        df = self.predict(pep_df)
        df["charge"] = self.charge_probs.apply(
            lambda x: self.charge_range[x>charge_prob]
        )
        df = df.explode("charge").dropna(subset=["charge"])
        if drop_probs_column:
            df.drop(columns="charge_probs", inplace=True)
        df["charge"] = df.charge.astype(np.int8)
        return df

class ChargeModelForAASeq(
    ModelInterface_for_Generic_AASeq_MultiLabelClassification
):
    """
    ModelInterface for charge prediction for amino acid sequence
    """
    def __init__(self, min_charge:int=1, max_charge:int=6):
        super().__init__(
            num_target_values=max_charge-min_charge+1,
            model_class=Model_for_Generic_AASeq_BinaryClassification_Transformer,
            nlayers=4, hidden_dim=128, dropout=0.1
        )

        self.target_column_to_predict = "charge_probs"
        self.target_column_to_train = "charge_indicators"
        self.min_charge = min_charge
        self.max_charge = max_charge
        self.charge_range = np.arange(min_charge, max_charge+1, dtype=np.int8)
        
    def predict_charges_for_pep_df(self, 
        pep_df:pd.DataFrame, 
        charge_prob=0.3,
        drop_probs_column=True
    ):
        df = self.predict(pep_df)
        df["charge"] = self.charge_probs.apply(
            lambda x: self.charge_range[x>charge_prob]
        )
        df = df.explode("charge").dropna(subset=["charge"])
        if drop_probs_column:
            df.drop(columns="charge_probs", inplace=True)
        df["charge"] = df.charge.astype(np.int8)
        return df

def group_psm_df_by_sequence(
    psm_df: pd.DataFrame,
    min_charge:int,
    max_charge:int,
):
    return psm_df.groupby("sequence")["charge"].apply(
        lambda x: get_charge_indicators(set(x),
            min_charge=min_charge, max_charge=max_charge
        )
    ).reset_index(drop=False).rename(columns={"charge":"charge_indicators"})


def group_psm_df_by_modseq(
    psm_df: pd.DataFrame,
    min_charge:int,
    max_charge:int,
):
    return psm_df.groupby(["sequence","mods","mod_sites"])["charge"].apply(
        lambda x: get_charge_indicators(set(x),
            min_charge=min_charge, max_charge=max_charge
        )
    ).reset_index(drop=False).rename(columns={"charge":"charge_indicators"})

def get_charge_indicators(
    charge_list,
    min_charge:int,
    max_charge:int,
):
    charge_indicators = np.zeros(max_charge-min_charge+1)
    for charge in charge_list:
        if charge <= max_charge and charge >= min_charge:
            charge_indicators[charge-min_charge] = 1.0
    return charge_indicators