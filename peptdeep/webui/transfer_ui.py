import streamlit as st
import pandas as pd
import os
import time

from peptdeep.settings import global_settings
from peptdeep.pipeline_api import transfer_learn
from peptdeep.webui.ui_utils import (
    files_in_pandas, select_files
)

def nce_search():
    grid_nce_first = st.number_input('Start NCE for grid NCE search',value = global_settings['model_mgr']['transfer']['grid_nce_first']*1.0,step = 1.0)
    global_settings['model_mgr']['transfer']['grid_nce_first'] = grid_nce_first
    grid_nce_last = st.number_input('End NCE for grid NCE search',min_value = grid_nce_first, value = global_settings['model_mgr']['transfer']['grid_nce_last']*1.0,step = 1.0)
    global_settings['model_mgr']['transfer']['grid_nce_last'] = grid_nce_last
    grid_nce_step = st.number_input('Step NCE for grid NCE search', value = global_settings['model_mgr']['transfer']['grid_nce_step']*1.0,step = 1.0)
    global_settings['model_mgr']['transfer']['grid_nce_step'] = grid_nce_step

    grid_instrument = st.multiselect('Instruments for grid NCE search', (
        global_settings['model_mgr']['instrument_group']
    ),default = global_settings['model_mgr']['transfer']['grid_instrument']) 
    global_settings['model_mgr']['transfer']['grid_instrument'] = grid_instrument

def fine_tune():
    epoch_ms2 = st.number_input('Epoch to train MS2 model', value = global_settings['model_mgr']['transfer']['epoch_ms2'])
    global_settings['model_mgr']['transfer']['epoch_ms2'] = epoch_ms2
    warmup_epoch_ms2 = st.number_input('Warmup epoch to train MS2 model', value = global_settings['model_mgr']['transfer']['epoch_ms2'], max_value=epoch_ms2)
    global_settings['model_mgr']['transfer']['warmup_epoch_ms2'] = warmup_epoch_ms2
    batch_size_ms2 = st.number_input('Mini-batch size to train MS2 model', value = global_settings['model_mgr']['transfer']['batch_size_ms2'])
    global_settings['model_mgr']['transfer']['batch_size_ms2'] = batch_size_ms2
    lr_ms2 = st.number_input('Learning rate to train MS2 model', value = global_settings['model_mgr']['transfer']['lr_ms2'], format='%e', step=1e-5)
    global_settings['model_mgr']['transfer']['lr_ms2'] = lr_ms2
    
    epoch_rt_ccs = st.number_input('Epoch to train RT and CCS models', value = global_settings['model_mgr']['transfer']['epoch_rt_ccs'])
    global_settings['model_mgr']['transfer']['epoch_rt_ccs'] = epoch_rt_ccs
    warmup_epoch_rt_ccs = st.number_input('Warmup epoch to train RT and CCS model', value = global_settings['model_mgr']['transfer']['epoch_rt_ccs'], max_value=epoch_rt_ccs)
    global_settings['model_mgr']['transfer']['warmup_epoch_rt_ccs'] = warmup_epoch_rt_ccs
    batch_size_rt_ccs = st.number_input('Mini-batch size to train RT and CCS model', value = global_settings['model_mgr']['transfer']['batch_size_rt_ccs'])
    global_settings['model_mgr']['transfer']['batch_size_rt_ccs'] = batch_size_rt_ccs
    lr_rt_ccs = st.number_input('Learning rate to train RT and CCS model', value = global_settings['model_mgr']['transfer']['lr_rt_ccs'], format='%e', step=1e-5)
    global_settings['model_mgr']['transfer']['lr_rt_ccs'] = lr_rt_ccs

def show():
    st.write("# Transfer model setup")

    model_output_folder = st.text_input('Model output folder')
    global_settings['model_mgr']['transfer']['model_output_folder'] = model_output_folder

    psm_type = st.selectbox('PSM type choice',global_settings['model_mgr']['transfer']['psm_type_choices'], index = 0)
    global_settings['model_mgr']['transfer']['psm_type'] = psm_type
    select_files(global_settings['model_mgr']['transfer']['psm_files'], "PSM files")
    ms_file_type = st.selectbox('MS file type',global_settings['model_mgr']['transfer']['ms_file_type_choices'], index = 0)
    global_settings['model_mgr']['transfer']['ms_file_type'] = ms_file_type
    ms_files = st.text_input('MS file folder')
    global_settings['model_mgr']['transfer']['ms_files'] = ms_files

    ms2_ppm = st.checkbox('MS2 ppm (otherwise Da)', global_settings['peak_matching']['ms2_ppm'])
    #ms2_ppm = st.selectbox('MS2 ppm',('True','False'))
    global_settings['peak_matching']['ms2_ppm'] = ms2_ppm
    ms2_tol_value = st.number_input('MS2 tolerance', value = global_settings['peak_matching']['ms2_tol_value'], step = 0.1)
    global_settings['peak_matching']['ms2_tol_value'] = ms2_tol_value

    fine_tune()

    psm_num_to_train_ms2 = st.number_input('PSM num to tune MS2 model', value = int(global_settings['model_mgr']['transfer']['psm_num_to_train_ms2']), step = 1)
    global_settings['model_mgr']['transfer']['psm_num_to_train_ms2'] = psm_num_to_train_ms2
    psm_num_per_mod_to_train_ms2 = st.number_input('PSM num per mod to tune MS2 model', value = int(global_settings['model_mgr']['transfer']['psm_num_per_mod_to_train_ms2']), step = 1)
    global_settings['model_mgr']['transfer']['psm_num_per_mod_to_train_ms2'] = psm_num_per_mod_to_train_ms2
    psm_num_to_train_rt_ccs = st.number_input('PSM num to tune RT and CCS model', value = int(global_settings['model_mgr']['transfer']['psm_num_to_train_rt_ccs']), step = 1)
    global_settings['model_mgr']['transfer']['psm_num_to_train_rt_ccs'] = psm_num_to_train_rt_ccs
    psm_num_per_mod_to_train_rt_ccs = st.number_input('PSM num per mod to tune RT and CCS model', value = int(global_settings['model_mgr']['transfer']['psm_num_per_mod_to_train_rt_ccs']), step = 1)
    global_settings['model_mgr']['transfer']['psm_num_per_mod_to_train_rt_ccs'] = psm_num_per_mod_to_train_rt_ccs
    top_n_mods_to_train = st.number_input('Top n mods to tune', value = int(global_settings['model_mgr']['transfer']['top_n_mods_to_train']), step = 1)
    global_settings['model_mgr']['transfer']['top_n_mods_to_train'] = top_n_mods_to_train


    st.write('### Grid NCE and instrument search')
    st.write('If NCE and instrument are unknown, grid search will look for the best values)')
    grid_nce_search = st.checkbox('Enabled', global_settings['model_mgr']['transfer']['grid_nce_search'])
    global_settings['model_mgr']['transfer']['grid_nce_search'] = grid_nce_search
    if grid_nce_search is True:
        nce_search()
    
    if st.button('Run transfer model'):
        transfer_learn()
        st.write('Transfer learning done!')