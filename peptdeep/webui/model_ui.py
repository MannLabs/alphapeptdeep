from email.policy import default
from faulthandler import disable
from multiprocessing.sharedctypes import Value
from tkinter import DISABLED
import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep
import pandas as pd
import os
import time
from peptdeep.settings import global_settings

def nce_search():
    grid_nce_first = st.number_input('Start NCE for grid NCE search',value = 15.0,step = 1.0)
    global_settings['model_mgr']['fine_tune']['grid_nce_first'] = grid_nce_first
    grid_nce_last = st.number_input('End NCE for grid NCE search',min_value = grid_nce_first, value = 45.0,step = 1.0)
    global_settings['model_mgr']['fine_tune']['grid_nce_last'] = grid_nce_last
    grid_nce_step = st.number_input('Step NCE for grid NCE search', value = 3.0,step = 1.0)
    global_settings['model_mgr']['fine_tune']['grid_nce_step'] = grid_nce_step

    grid_instrument = st.multiselect('Instruments for grid NCE search', (
        global_settings['model_mgr']['predict']['instrument_choices']
    ),default =['Lumos']) 
    global_settings['model_mgr']['fine_tune']['grid_instrument'] = grid_instrument

def fine_tune():
    epoch_ms2 = st.number_input('Epoch to tune MS2 model', value = 10)
    global_settings['model_mgr']['fine_tune']['epoch_ms2'] = epoch_ms2
    epoch_rt_ccs = st.number_input('Epoch to tune RT and CCS models', value = 20)
    global_settings['model_mgr']['fine_tune']['epoch_rt_ccs'] = epoch_rt_ccs


def predict():
    batch_size_ms2 = st.number_input('Batch size to tune MS2 model', value = 512)
    global_settings['model_mgr']['predict']['batch_size_ms2'] = batch_size_ms2
    batch_size_rt_ccs = st.number_input('Batch size to tune RT and CCS model', value = 1024)
    global_settings['model_mgr']['predict']['batch_size_rt_ccs'] = batch_size_rt_ccs
    default_nce = st.number_input('NCE', value = 30)
    global_settings['model_mgr']['predict']['default_nce'] = default_nce
    default_instrument = st.selectbox('Instrument',(global_settings['model_mgr']['predict']['instrument_choices']),index = 0)
    global_settings['model_mgr']['predict']['default_instrument'] = default_instrument

    verbose = st.checkbox('Verbose')
    global_settings['model_mgr']['predict']['verbose'] = verbose
    multiprocessing = st.checkbox('Multiprocessing')
    global_settings['model_mgr']['predict']['multiprocessing'] = multiprocessing

def model():
    model_url = st.text_input('URL to download the pre-trained models:',value = 'https://datashare.biochem.mpg.de/s/ABnQuD2KkXfIGF3/download')
    global_settings['model_url'] = model_url

    model_url_zip_name = st.text_input('The downloaded zip file name from the URL:',value = 'released_models')
    global_settings['model_url_zip_name'] = model_url_zip_name

    thread_num = st.number_input('thread number:', value = 8)
    global_settings['thread_num'] = thread_num


def show():

    st.write("# Model Configuration")
    st.write('### Model parameters')
    model()

    model_type = st.selectbox('Model type',(global_settings['model_mgr']['model_choices']),index = 0)
    global_settings['model_mgr']['model_type'] = model_type

    st.write('### Prediction parameters')
    predict()

    st.write('### Fine-tuning parameters')
    fine_tune()

    grid_nce_search = st.checkbox('grid NCE search')
    global_settings['model_mgr']['fine_tune']['grid_nce_search'] = grid_nce_search
    if grid_nce_search is True:
        nce_search()

    global_settings['model_mgr']['mask_modloss'] = bool(st.checkbox('mask modloss',value = True))
