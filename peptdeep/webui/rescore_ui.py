import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep
import pandas as pd
import os
import time

def files_in_folder_pandas(folder: str) -> pd.DataFrame:
    """Reads a folder and returns a pandas dataframe containing the files and additional information.
    Args:
        folder (str): Path to folder.

    Returns:
        pd.DataFrame: PandasDataFrame.
    """
    files = os.listdir(folder)
    created = [time.ctime(os.path.getctime(os.path.join(folder, _))) for _ in files]
    sizes = [os.path.getsize(os.path.join(folder, _)) / 1024 ** 2 for _ in files]
    df = pd.DataFrame(files, columns=["File"])
    df["Created"] = created
    df["Filesize (Mb)"] = sizes

    return df

def show():
    """Streamlit page that displays information on how to get started."""
    st.write("# Rescore started")
    st.text("Welcome to rescore of PeptDeep.")

    markdown_link("google link", "https://www.google.com")
    fasta_files = files_in_folder_pandas("/Users/zhouxiexuan/workspace/alphadeep/alphadeep")

    st.table(fasta_files)
    option_PSM_type = st.selectbox(
     'PSM types',
     ('pFind', 'AlphaPept', 'MaxQuant'))
    #st.write('You selected:', option_PSM_type)
    option_model_type = st.selectbox(
        'Model type',
        ('regular','phos','HLA')
    )
    agree = st.checkbox('I agree')
    if agree:
     st.write('Great!')