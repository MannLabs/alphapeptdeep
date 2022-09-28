import yaml
import streamlit as st
from io import StringIO
import multiprocessing

from alphabase.constants.modification import (
    MOD_DF, keep_modloss_by_importance
)

from peptdeep.settings import (
    global_settings,
    update_settings,
)

def show():
    load_settings_gui()
    save_settings_gui()

    st.write("### Common settings")

    ms2_ppm = st.checkbox(label='MS2 ppm (otherwise Da)', value=global_settings['peak_matching']['ms2_ppm'])
    global_settings['peak_matching']['ms2_ppm'] = ms2_ppm
    ms2_tol_value = st.number_input(label='MS2 tolerance', value = global_settings['peak_matching']['ms2_tol_value'], step = 1.0)
    global_settings['peak_matching']['ms2_tol_value'] = ms2_tol_value

    ms1_ppm = st.checkbox(label='MS1 ppm (otherwise Da)', value=global_settings['peak_matching']['ms1_ppm'])
    global_settings['peak_matching']['ms1_ppm'] = ms1_ppm
    ms1_tol_value = st.number_input(label='MS1 tolerance', value = global_settings['peak_matching']['ms1_tol_value'], step = 1.0)
    global_settings['peak_matching']['ms1_tol_value'] = ms1_tol_value
    
    thread_num = st.number_input(label='Thread number', 
        value=global_settings['thread_num'], 
        max_value=multiprocessing.cpu_count(), step=1
    )
    global_settings['thread_num'] = thread_num

    device_type = st.selectbox(label='Computing devices',
        options=global_settings['torch_device']['device_type_choices'],
        index = global_settings['torch_device']['device_type_choices'].index(
            global_settings['torch_device']['device_type']
        )
    )
    global_settings['torch_device']['device_type'] = device_type

    log_level = st.selectbox(label='Log level', 
        options=global_settings['log_level_choices'], 
        index = global_settings['log_level_choices'].index(
            global_settings['log_level']
        )
    )
    global_settings['log_level'] = log_level

    global_settings['common']['modloss_importance_level'] = st.number_input(
        'Modification loss importance level (for a PTM, fragment modloss mz=0 if modloss_importance<modloss_importance_level)', 
        value=global_settings['common']['modloss_importance_level'], step=1.0,
    )
    keep_modloss_by_importance(global_settings['common']['modloss_importance_level'])

    st.write("Modification modloss example (check the `modloss` column):")
    st.dataframe(MOD_DF.loc[
        ['Carbamidomethyl@C','Oxidation@M','Phospho@S'],
        ['mod_name','classification','mass','modloss','modloss_original','modloss_importance']
    ])

def _update_st_session_state_after_loading_settings(
    state_dict:dict={
        'select_psm_type': global_settings['model_mgr']['transfer']['psm_type'],
        'select_ms_file_type': global_settings['model_mgr']['transfer']['ms_file_type'],
        'lib_input_type': global_settings['library']['input']['infile_type'],
    }
):
    for key, val in state_dict.items():
        st.session_state[key] = val

def load_settings_gui():
    st.write("### Load previous yaml settings")
    uploaded_file = st.file_uploader("Choose a local yaml file")
    if uploaded_file is not None:
        f = StringIO(uploaded_file.getvalue().decode("utf-8"))
        uploaded_settings = yaml.load(f, Loader=yaml.FullLoader)
        update_settings(global_settings, uploaded_settings)
        st.write("Global settings have been updated")

def save_settings_gui():
    st.write("### Save current settings")
    st.write("The saved yaml file can be used as a template for CLI commands)")

    f = StringIO()
    yaml.dump(global_settings, f)

    st.download_button(
        label="Download settings as yaml",
        data=f.getvalue(),
        file_name='peptdeep_settings.yaml',
    )