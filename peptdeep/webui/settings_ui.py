import yaml
import streamlit as st
import pandas as pd
from io import StringIO
import multiprocessing

from alphabase.constants.modification import (
    MOD_DF, keep_modloss_by_importance
)

from peptdeep.settings import (
    global_settings,
    update_settings,
    add_user_defined_modifications,
)

def add_user_mods():
    st.write("#### User-defined modifications")
    st.write('PeptDeep supports modifications those are not in UniMod')
    user_mod_expander = st.expander(label="Add user-defined modifications")
    with user_mod_expander:
        mod_name = st.text_input(
            label='User-defined modification name, e.g. Hello@K',
            key='user_mod_name'
        ).strip()
        composition = st.text_input(
            label='The modification composition, e.g. H(1)P(1)O(3)',
            key='user_mod_comp'
        ).strip()
        modloss_composition = st.text_input(
            label="The modification loss composition, e.g. H(3)P(1)O(4)",
            key='user_mod_loss'
        ).strip()

        if mod_name:
            global_settings['common']['user_defined_modifications'][mod_name] = {
                'composition': composition,
                'modloss_composition': modloss_composition,
            }

        st.dataframe(pd.DataFrame().from_dict(
            global_settings['common']['user_defined_modifications'],
            orient = 'index',
        ))

        def _clear_user_mods():
            global_settings['common']['user_defined_modifications'] = {}
            st.session_state.user_mod_name = ''
            st.session_state.user_mod_comp = ''
            st.session_state.user_mod_loss = ''

        st.button(label='Clear all user modifications', 
            on_click=_clear_user_mods
        )

        if st.button(label='Add user modifications into AlphaBase'):
            add_user_defined_modifications()
            st.write("Check last n+2 modifications:")
            st.dataframe(MOD_DF.tail(
                len(global_settings['common'][
                    'user_defined_modifications'
                ])+2
            ))

def show():
    load_settings_gui()
    save_settings_gui()

    st.write("### Common settings")
    
    add_user_mods()

    ms2_ppm = st.checkbox(label='MS2 ppm (otherwise Da)', value=global_settings['peak_matching']['ms2_ppm'])
    global_settings['peak_matching']['ms2_ppm'] = ms2_ppm
    ms2_tol_value = st.number_input(label='MS2 tolerance', value = global_settings['peak_matching']['ms2_tol_value'], step = 1.0)
    global_settings['peak_matching']['ms2_tol_value'] = ms2_tol_value

    ms1_ppm = st.checkbox(label='MS1 ppm (otherwise Da)', value=global_settings['peak_matching']['ms1_ppm'])
    global_settings['peak_matching']['ms1_ppm'] = ms1_ppm
    ms1_tol_value = st.number_input(label='MS1 tolerance', value = global_settings['peak_matching']['ms1_tol_value'], step = 1.0)
    global_settings['peak_matching']['ms1_tol_value'] = ms1_tol_value
    
    cpu_count = multiprocessing.cpu_count()
    thread_num = st.number_input(label='Thread number', 
        value=global_settings['thread_num'] 
          if global_settings['thread_num']<cpu_count 
          else cpu_count, 
        max_value=cpu_count, step=1
    )
    global_settings['thread_num'] = thread_num

    global_settings['torch_device']['device_type'] = st.selectbox(
        label='Computing devices',
        options=global_settings['torch_device']['device_type_choices'],
        index = global_settings['torch_device']['device_type_choices'].index(
            global_settings['torch_device']['device_type']
        )
    )

    global_settings['log_level'] = st.selectbox(
        label='Log level', 
        options=global_settings['log_level_choices'], 
        index = global_settings['log_level_choices'].index(
            global_settings['log_level']
        )
    )

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
        'lib_input_type': global_settings['library']['infile_type'],
    }
):
    for key, val in state_dict.items():
        st.session_state[key] = val

def load_settings_gui():
    st.write("### Load previous yaml settings")
    uploaded_file = st.file_uploader("Choose a yaml file")
    if uploaded_file is not None:
        f = StringIO(uploaded_file.getvalue().decode("utf-8"))
        uploaded_settings = yaml.load(f, Loader=yaml.FullLoader)
        update_settings(global_settings, uploaded_settings)
        st.write("Global settings have been updated")

def save_settings_gui():
    st.write("### Save current settings")
    st.write("The saved yaml file can be used in command line interface")

    f = StringIO()
    yaml.dump(global_settings, f)

    st.download_button(
        label="Download settings as yaml",
        data=f.getvalue(),
        file_name='peptdeep_settings.yaml',
    )