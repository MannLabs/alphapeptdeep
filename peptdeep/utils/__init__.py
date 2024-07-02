from .logger import *
from .regression import *
from .device_utils import *

import os
import tqdm
import logging

import itertools

import pandas as pd
import numpy as np


# from alphatims
def process_bar(iterator, len_iter):
    with tqdm.tqdm(total=len_iter) as bar:
        i = 0
        for i, iter in enumerate(iterator):
            yield iter
            bar.update()
        bar.update(len_iter - i - 1)


def _get_delimiter(tsv_file: str):
    with open(tsv_file, "r") as f:
        line = f.readline().strip()
        if "\t" in line:
            return "\t"
        elif "," in line:
            return ","
        else:
            return "\t"


def read_peptide_table(tsv_file: str) -> pd.DataFrame:
    sep = _get_delimiter(tsv_file)
    df = pd.read_csv(tsv_file, sep=sep, keep_default_na=False)
    if "mod_sites" in df.columns:
        df["mod_sites"] = df.mod_sites.astype("U")
    return df


_special_raw_suffices = [".ms_data.hdf", "_hcdft.mgf", ".mzml" ".mgf"]


def parse_ms_file_names_to_dict(ms_file_list: list) -> dict:
    """
    Load spectrum file paths into a dict:
        "/Users/xxx/raw_name.raw" -> {"raw_name":"/Users/xxx/raw_name.raw"}

    Parameters
    ----------
    spectrum_file_list : list
        File path list

    Returns
    -------
    dict
        {"raw_name":"/Users/xxx/raw_name.raw", ...}
    """

    spec_dict = {}
    for ms_file in ms_file_list:
        raw_filename = os.path.basename(ms_file)
        raw_name = raw_filename.lower()
        for raw_suff in _special_raw_suffices:
            if raw_name.endswith(raw_suff):
                raw_name = raw_filename[: -len(raw_suff)]
                break
        if len(raw_filename) == len(raw_name):
            raw_name = os.path.splitext(raw_name)[0]
        spec_dict[raw_name] = ms_file
    return spec_dict


def _flatten(list_of_lists):
    """
    Flatten a list of lists
    """
    return list(itertools.chain.from_iterable(list_of_lists))


def explode_multiple_columns(df: pd.DataFrame, columns: list):
    try:
        return df.explode(columns)
    except ValueError:
        # pandas <= 1.2.x?
        logging.warn(f"pandas=={pd.__version__} cannot explode multiple columns")
        ret_df = df.explode(columns[0])
        for col in columns[1:]:
            ret_df[col] = _flatten(df[col].values)
        return ret_df
