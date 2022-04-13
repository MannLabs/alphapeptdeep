from email.policy import default
from multiprocessing.sharedctypes import Value
import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep
import pandas as pd
import os
import time
from peptdeep.settings import global_settings

def show():

    st.write("# Model ")
    st.text("Welcome to model of AlphaPeptDeep.")
    st.write('### model')
    model_url = st.text_input('model_url:',value = 'https://datashare.biochem.mpg.de/s/ABnQuD2KkXfIGF3/download')
    global_settings['model_url'] = model_url

    model_url_zip_name = st.text_input('model url zip name:',value = 'released_models')
    global_settings['model_url_zip_name'] = model_url_zip_name

    local_model_file_name = st.text_input('local model zip name:',value = 'pretrained_models.zip')
    global_settings['local_model_zip_name'] = local_model_file_name

    thread_num = st.number_input('thread num:', value = 8)
    global_settings['thread_num'] = thread_num

    st.write('predict')
    batch_size_ms2 = st.number_input('batch size ms2', value = 512)
    global_settings['model_mgr']['predict']['batch_size_ms2'] = batch_size_ms2
    batch_size_rt_ccs = st.number_input('batch size rt ccs', value = 1024)
    global_settings['model_mgr']['predict']['batch_size_rt_ccs'] = batch_size_rt_ccs
    default_nce = st.number_input('nce', value = 30)
    global_settings['model_mgr']['predict']['default_nce'] = default_nce
    default_instrument = st.selectbox('instrument',(global_settings['model_mgr']['predict']['instrument_choices']),index = 0)
    global_settings['model_mgr']['predict']['default_instrument'] = default_instrument

    verbose = st.selectbox('verbose',('True','False'),index = 0)
    global_settings['model_mgr']['predict']['verbose'] = verbose
    multiprocessing = st.selectbox('multiprocessing',('True','False'),index = 0)
    global_settings['model_mgr']['predict']['multiprocessing'] = multiprocessing


    epoch_ms2 = st.number_input('epoch ms2', value = 10)
    global_settings['model_mgr']['fine_tune']['epoch_ms2'] = epoch_ms2
    epoch_rt_ccs = st.number_input('epoch rt ccs', value = 20)
    global_settings['model_mgr']['fine_tune']['epoch_rt_ccs'] = epoch_rt_ccs
    grid_nce_search = st.selectbox('grid nce search',('True','False'),index = 0)
    global_settings['model_mgr']['fine_tune']['grid_nce_search'] = grid_nce_search

    grid_nce_first = st.number_input('grid nce first',value = 15.0,step = 1.0)
    global_settings['model_mgr']['fine_tune']['grid_nce_first'] = grid_nce_first
    grid_nce_last = st.number_input('grid nce last',value = 45.0,step = 1.0)
    global_settings['model_mgr']['fine_tune']['grid_nce_last'] = grid_nce_last
    grid_nce_step = st.number_input('grid nce step', value = 3.0,step = 1.0)
    global_settings['model_mgr']['fine_tune']['grid_nce_step'] = grid_nce_step

    grid_instrument = st.multiselect('grid instrument', (
        global_settings['model_mgr']['predict']['instrument_choices']
    ),default =['Lumos']) 
    global_settings['model_mgr']['fine_tune']['grid_instrument'] = grid_instrument


    global_settings['model_mgr']['mask_modloss'] = bool(st.checkbox('mask modloss'))
    model_type = st.selectbox('model choices',(global_settings['model_mgr']['model_choices']),index = 0)
    global_settings['model_mgr']['model_type'] = model_type
