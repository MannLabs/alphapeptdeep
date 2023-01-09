import os
import collections

from alphabase.yaml_utils import load_yaml
from alphabase.constants.modification import (
    load_mod_df, keep_modloss_by_importance,
    add_new_modifications
)

from peptdeep.constants._const import CONST_FOLDER

global_settings = load_yaml(
    os.path.join(CONST_FOLDER, 'default_settings.yaml')
)
"""
Global settings in peptdeep, 
it controls all functionalities of PeptDeep.
"""

def _refine_global_settings():
    global_settings['PEPTDEEP_HOME'] = os.path.expanduser(
        global_settings['PEPTDEEP_HOME']
    )
    global_settings['library']['output_folder']=(
        global_settings['library']['output_folder'].format(
            PEPTDEEP_HOME=global_settings['PEPTDEEP_HOME']
        )
    )
    global_settings['model_mgr']['transfer']['model_output_folder']=(
        global_settings['model_mgr']['transfer']['model_output_folder'].format(
            PEPTDEEP_HOME=global_settings['PEPTDEEP_HOME']
        )
    )
    global_settings['percolator']['output_folder']=(
        global_settings['percolator']['output_folder'].format(
            PEPTDEEP_HOME=global_settings['PEPTDEEP_HOME']
        )
    )
    for key, val in list(global_settings['model_mgr'][
        'instrument_group'
    ].items()):
        global_settings['model_mgr'][
            'instrument_group'
        ][key.upper()] = val
_refine_global_settings()

model_const = load_yaml(
    os.path.join(CONST_FOLDER, 'model_const.yaml')
)

def update_settings(dict_, new_dict):
    for k, v in new_dict.items():
        if isinstance(v, collections.abc.Mapping):
            dict_[k] = update_settings(dict_.get(k, {}), v)
        else:
            dict_[k] = v
    return dict_

def load_global_settings(yaml:str):
    d = load_yaml(yaml)
    update_settings(global_settings, d)
    _refine_global_settings()

def update_modifications(tsv:str="", 
    modloss_importance_level:float=1.0
):
    """
    Load modification tsv either from alphabase default 
    `modification.tsv <https://github.com/MannLabs/alphabase/blob/main/alphabase/constants/const_files/modification.tsv>`_
    or an external tsv file.

    Parameters
    ----------
    tsv : str, optional
        External modification tsv file, "" refers to the default alphabase tsv,
        by default "".
    modloss_importance_level : float, optional
        Only keep the important modification losses, by default 1.0
    """
    if os.path.isfile(tsv):
        load_mod_df(tsv, modloss_importance_level=modloss_importance_level)
    else:
        keep_modloss_by_importance(modloss_importance_level)
    
    add_user_defined_modifications()

def add_user_defined_modifications(
    user_mods:dict=None
):
    """
    Add user-defined modifications into the system,
    this is userful for isotope labeling. 

    Parameters
    ----------
    user_mods : dict, optional
        Example:
        ```
        {
        "Dimethyl2@Any N-term": { 
        "composition": "H(2)2H(2)C(2)",
        "modloss_composition": ""
        }, ...
        }
        ```
        Set as `global_settings["user_defined_modifications"]` if it is None.
        By default None.
    """
    if user_mods is None:
        user_mods = global_settings["common"]["user_defined_modifications"]
    add_new_modifications(user_mods)

    from peptdeep.model.featurize import update_all_mod_features
    update_all_mod_features()

update_modifications()
