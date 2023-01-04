import streamlit as st
import pandas as pd
import os
import time

from peptdeep.settings import global_settings
from peptdeep.webui.ui_utils import files_in_folder_pandas

def show():
    """Streamlit page that displays information on how to rescore."""
    st.write("# DDA Rescoring")

    MS_type = st.selectbox(
        label='MS file type',
        options=('Raw', 'MGF', 'ms_data.hdf')
    )
    raw_folder = st.text_input(label='Raw folder')
    if os.path.isdir(raw_folder):
        st.write(f"### MS files in {raw_folder}")

        raw_files = files_in_folder_pandas(raw_folder,MS_type)

        st.dataframe(raw_files)
    else:
        st.write(f"Invalid folder: {raw_folder}")

    result_folder = st.text_input(label='Result folder')
    #st.write('The current result folder is', result_folder)
    PSM_type = st.selectbox(label='PSM file type',
     options=('AlphaPept', 'pFind')
    )
    if PSM_type == 'AlphaPept':
        psm_type = 'ms_data.hdf'
    elif PSM_type == 'pFind':
        psm_type = 'spectra'

    if os.path.isdir(result_folder):

        st.write(f"### PSM files in {result_folder}")

        result_files = files_in_folder_pandas(result_folder,psm_type)

        st.dataframe(result_files)

    st.warning("We are still working on Rescore GUI panel.")
    st.warning("For command line users, please use `peptdeep rescore` for rescoring.")
    st.warning("For Python users, please use `peptdeep.pipeline_api.rescore()` for rescoring.")