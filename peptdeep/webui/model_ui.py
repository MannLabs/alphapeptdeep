import streamlit as st
from peptdeep.settings import global_settings
import multiprocessing

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
    epoch_ms2 = st.number_input('Epoch to tune MS2 model', value = global_settings['model_mgr']['transfer']['epoch_ms2'])
    global_settings['model_mgr']['transfer']['epoch_ms2'] = epoch_ms2
    epoch_rt_ccs = st.number_input('Epoch to tune RT and CCS models', value = global_settings['model_mgr']['transfer']['epoch_rt_ccs'])
    global_settings['model_mgr']['transfer']['epoch_rt_ccs'] = epoch_rt_ccs


def predict():
    batch_size_ms2 = st.number_input('Batch size to predict MS2', value = global_settings['model_mgr']['predict']['batch_size_ms2'])
    global_settings['model_mgr']['predict']['batch_size_ms2'] = batch_size_ms2
    batch_size_rt_ccs = st.number_input('Batch size to predict RT and CCS', value = global_settings['model_mgr']['predict']['batch_size_rt_ccs'])
    global_settings['model_mgr']['predict']['batch_size_rt_ccs'] = batch_size_rt_ccs

    default_instrument = st.selectbox('Instrument',(list(global_settings['model_mgr']['instrument_group'].keys())),index = 0)
    global_settings['model_mgr']['default_instrument'] = default_instrument
    default_nce = st.number_input('NCE', value = global_settings['model_mgr']['default_nce'],disabled=(default_instrument=='timsTOF'))
    global_settings['model_mgr']['default_nce'] = default_nce


    verbose = st.checkbox('Verbose', global_settings['model_mgr']['predict']['verbose'])
    global_settings['model_mgr']['predict']['verbose'] = verbose
    multiprocessing = st.checkbox('Multiprocessing', global_settings['model_mgr']['predict']['multiprocessing'])
    global_settings['model_mgr']['predict']['multiprocessing'] = multiprocessing

def model():
    model_url = st.text_input('URL (or local path) to download the pre-trained models',value = global_settings['model_url'])
    global_settings['model_url'] = model_url
    
    thread_num = st.number_input('Thread number', value = multiprocessing.cpu_count()-1)
    global_settings['thread_num'] = thread_num

    global_settings['model_mgr']['external_ms2_model'] = st.text_input('External MS2 model')
    global_settings['model_mgr']['external_rt_model'] = st.text_input('External RT model')
    global_settings['model_mgr']['external_ccs_model'] = st.text_input('External CCS model')


def show():

    st.write("# Model Configuration")
    st.write('### Model parameters')
    model()

    model_type = st.selectbox('Model type',(global_settings['model_mgr']['model_choices']),index = 0)
    global_settings['model_mgr']['model_type'] = model_type
    global_settings['model_mgr']['mask_modloss'] = bool(
        st.checkbox('mask modloss (setting intensity values to zero for neutral loss of PTMs (e.g. -98 Da for Phospho@S/T))',
        value = global_settings['model_mgr']['mask_modloss'])
    )

    st.write('### Prediction parameters')
    predict()

    st.write('### Fine-tuning parameters (for DDA rescoring only)')
    fine_tune()

    st.write('### Grid NCE and instrument search for DDA rescoring')
    grid_nce_search = st.checkbox('Enabled', global_settings['model_mgr']['transfer']['grid_nce_search'])
    global_settings['model_mgr']['transfer']['grid_nce_search'] = grid_nce_search
    if grid_nce_search is True:
        nce_search()
