import streamlit as st
from peptdeep.settings import global_settings

def predict():
    batch_size_ms2 = st.number_input(label='Batch size to predict MS2', value = global_settings['model_mgr']['predict']['batch_size_ms2'])
    global_settings['model_mgr']['predict']['batch_size_ms2'] = batch_size_ms2
    batch_size_rt_ccs = st.number_input(label='Batch size to predict RT and CCS', value = global_settings['model_mgr']['predict']['batch_size_rt_ccs'])
    global_settings['model_mgr']['predict']['batch_size_rt_ccs'] = batch_size_rt_ccs

    instruments = list(global_settings['model_mgr']['instrument_group'].keys())
    global_settings['model_mgr']['default_instrument'] = st.selectbox(
        label='Instrument',options=instruments,index = instruments.index(
            global_settings['model_mgr']['default_instrument']
        )
    )
    default_nce = st.number_input(label='NCE', value = global_settings['model_mgr']['default_nce'])
    global_settings['model_mgr']['default_nce'] = default_nce

    verbose = st.checkbox(label='Verbose', value=global_settings['model_mgr']['predict']['verbose'])
    global_settings['model_mgr']['predict']['verbose'] = verbose
    multiprocessing = st.checkbox(label='Multiprocessing (if no GPUs)', value=global_settings['model_mgr']['predict']['multiprocessing'])
    global_settings['model_mgr']['predict']['multiprocessing'] = multiprocessing

def model():
    model_url = st.text_input(label='URL (or local path) to download the pre-trained models',value = global_settings['model_url'])
    global_settings['model_url'] = model_url

    global_settings['model_mgr']['external_ms2_model'] = st.text_input(label='External MS2 model', value=global_settings['model_mgr']['external_ms2_model'])
    global_settings['model_mgr']['external_rt_model'] = st.text_input(label='External RT model', value=global_settings['model_mgr']['external_rt_model'])
    global_settings['model_mgr']['external_ccs_model'] = st.text_input(label='External CCS model', value=global_settings['model_mgr']['external_ccs_model'])


def show():

    st.write("# Model Configuration")
    st.write('### Pre-trained models')
    model()

    global_settings['model_mgr']['model_type'] = st.selectbox(
        label='Model type',
        options=global_settings['model_mgr']['model_choices'],
        index = global_settings['model_mgr']['model_choices'].index(
            global_settings['model_mgr']['model_type']
        )
    )
    global_settings['model_mgr']['mask_modloss'] = bool(
        st.checkbox(label='mask modloss (this will set intensity values to zero for neutral loss of PTMs (e.g. -98 Da for Phospho@S/T))',
        value = global_settings['model_mgr']['mask_modloss'])
    )

    st.write('### Prediction')
    predict()
