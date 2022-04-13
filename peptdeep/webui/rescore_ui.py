import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep
import pandas as pd
import os
import time

from peptdeep.settings import global_settings

def files_in_folder_pandas(folder: str, file_type:str=None) -> pd.DataFrame:
    """Reads a folder and returns a pandas dataframe containing the files and additional information.
    Args:
        folder (str): Path to folder.

    Returns:
        pd.DataFrame: PandasDataFrame.
    """
    if file_type is None:
        files = os.listdir(folder)
    else:
        file_type = file_type.lower()
        files = [
            file for file in os.listdir(folder) 
            if file.lower().endswith(f".{file_type}")
        ]
    created = [time.ctime(os.path.getctime(os.path.join(folder, _))) for _ in files]
    sizes = [os.path.getsize(os.path.join(folder, _)) / 1024 ** 2 for _ in files]
    df = pd.DataFrame(files, columns=["File"])
    df["Created"] = created
    df["Filesize (Mb)"] = sizes

    return df

def show():
    """Streamlit page that displays information on how to rescore."""
    st.write("# Rescore started")
    st.text("Welcome to rescore of AlphaPeptDeep.")

    raw_folder = st.text_input('Raw folder')
    #st.write('The current raw folder is', raw_folder)
    MS_type = st.selectbox(
     'MS types',
     ('Raw', 'MGF', 'hdf','py'))
    #st.write('You selected:', MS_type)
    if raw_folder:
        st.text(
            f"PeptDeep looks for MS files in {raw_folder}.\nThese can be selected in the new experiment tab.\nYou can add own files to this folder."
            )

        st.write("### Existing files")

        raw_files = files_in_folder_pandas(raw_folder,MS_type)

        st.table(raw_files)

    result_folder = st.text_input('Result folder')
    #st.write('The current result folder is', result_folder)
    PSM_type = st.selectbox(
     'PSM types',
     ('AlphaPept', 'pFind', 'MaxQuant'))
    if PSM_type == 'AlphaPept':
        psm_type = 'hdf'
    if PSM_type == 'pFind':
        psm_type = 'spectra'
    if PSM_type == 'MaxQuant':
        psm_type = 'txt'
    #st.write('You selected:', PSM_type)
    if result_folder:
        st.text(
            f"PeptDeep looks for PSM files in {result_folder}.\nThese can be selected in the new experiment tab.\nYou can add own files to this folder."
            )

        st.write("### Existing files")

        result_files = files_in_folder_pandas(result_folder,psm_type)

        st.table(result_files)