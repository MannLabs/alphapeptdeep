import streamlit as st
import pandas as pd
import os
import time

from peptdeep.settings import global_settings

def show():
    st.write("# Transfer model")

    ms2_ppm = st.checkbox('MS2 ppm', global_settings['peak_matching']['ms2_ppm'])
    global_settings['peak_matching']['ms2_ppm'] = ms2_ppm
    ms2_tol_value = st.number_input('MS2 tol value', value = global_settings['peak_matching']['ms2_tol_value'], step = 0.5)
    global_settings['peak_matching']['ms2_tol_value'] = ms2_tol_value

    epoch_ms2 = st.number_input('Epoch to fine-tune MS2', value = global_settings['model_mgr']['fine_tune']['epoch_ms2'], step = 1)
    global_settings['model_mgr']['fine_tune']['epoch_ms2'] = epoch_ms2
    epoch_rt_ccs = st.number_input('Epoch RT CCS', value = global_settings['model_mgr']['fine_tune']['epoch_rt_ccs'], step = 1)
    global_settings['model_mgr']['fine_tune']['epoch_rt_ccs'] = epoch_rt_ccs
    grid_nce_search = st.checkbox('Grid nce search', global_settings['model_mgr']['fine_tune']['grid_nce_search'])
    global_settings['model_mgr']['fine_tune']['grid_nce_search'] = grid_nce_search

    grid_nce_first = st.number_input('Grid nce first', value = global_settings['model_mgr']['fine_tune']['grid_nce_first'], step = 1)
    global_settings['model_mgr']['fine_tune']['grid_nce_first'] = grid_nce_first
    grid_nce_last = st.number_input('Grid nce last', value = global_settings['model_mgr']['fine_tune']['grid_nce_last'], step = 1)
    global_settings['model_mgr']['fine_tune']['grid_nce_last'] = grid_nce_last
    grid_nce_step = st.number_input('Grid nce steps', value = global_settings['model_mgr']['fine_tune']['grid_nce_step'], step = 1)
    global_settings['model_mgr']['fine_tune']['grid_nce_step'] = grid_nce_step
    grid_instrument = st.selectbox('Grid instrument',global_settings['model_mgr']['instrument_group'])
    global_settings['model_mgr']['fine_tune']['grid_instrument'] = grid_instrument

    ms2_output_path = st.text_input('MS2 output path')
    global_settings['model_mgr']['fine_tune']['ms2_output_path'] = ms2_output_path
    rt_output_path = st.text_input('RT output path')
    global_settings['model_mgr']['fine_tune']['rt_output_path'] = rt_output_path
    ccs_output_path = st.text_input('CCS output path')
    global_settings['model_mgr']['fine_tune']['ccs_output_path'] = ccs_output_path

    psm_type = st.selectbox('PSM type choices',global_settings['model_mgr']['fine_tune']['psm_type_choices'], index = 0)
    global_settings['model_mgr']['fine_tune']['psm_type'] = psm_type
    psm_files = st.text_input('PSM file folder')
    global_settings['model_mgr']['fine_tune']['psm_files'] = psm_files
    ms_file_type = st.selectbox('MS file type',global_settings['model_mgr']['fine_tune']['ms_file_type_choices'], index = 0)
    global_settings['model_mgr']['fine_tune']['ms_file_type'] = ms_file_type
    ms_files = st.text_input('MS file folder')
    global_settings['model_mgr']['fine_tune']['ms_files'] = ms_files

    psm_num_to_tune_ms2 = st.number_input('PSM num to tune MS2', global_settings['model_mgr']['fine_tune']['psm_num_to_tune_ms2'], step = 1)
    global_settings['model_mgr']['fine_tune']['psm_num_to_tune_ms2'] = psm_num_to_tune_ms2
    psm_num_per_mod_to_tune_ms2 = st.number_input('PSM num per mod to tune MS2', global_settings['model_mgr']['fine_tune']['psm_num_per_mod_to_tune_ms2'], step = 1)
    global_settings['model_mgr']['fine_tune']['psm_num_per_mod_to_tune_ms2'] = psm_num_per_mod_to_tune_ms2
    psm_num_to_tune_rt_ccs = st.number_input('PSM num to tune RT CCS', global_settings['model_mgr']['fine_tune']['psm_num_to_tune_rt_ccs'], step = 1)
    global_settings['model_mgr']['fine_tune']['psm_num_to_tune_rt_ccs'] = psm_num_to_tune_rt_ccs
    psm_num_per_mod_to_tune_rt_ccs = st.number_input('PSM num per mod to tune RT CCS', global_settings['model_mgr']['fine_tune']['psm_num_per_mod_to_tune_rt_ccs'], step = 1)
    global_settings['model_mgr']['fine_tune']['psm_num_per_mod_to_tune_rt_ccs'] = psm_num_per_mod_to_tune_rt_ccs
    top_n_mods_to_tune = st.number_input('Top n mods to tune', global_settings['model_mgr']['fine_tune']['top_n_mods_to_tune'], step = 1)
    global_settings['model_mgr']['fine_tune']['top_n_mods_to_tune'] = top_n_mods_to_tune