import streamlit as st
import pandas as pd
import os
import time

from peptdeep.settings import global_settings

def show():
    st.write("# Transfer model")
    ms2_output_path = st.text_input('ms2_output_path')
    rt_output_path = st.text_input('rt_output_path')
    ccs_output_path = st.text_input('ccs_output_path')
    psm_type = st.selectbox('psm_type:',('alphapept','pfind','maxquant','diann'))
    
    psm_files = st.text_input('psm_files')
    ms_file_type = st.selectbox('ms_file_type',('alphapept_hdf','thermo_raw','mgf','mzml'))
    
    ms_files = st.text_input('ms_files')