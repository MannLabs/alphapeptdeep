import streamlit as st
import pandas as pd
import os
import time
from datetime import datetime

from alphabase.yaml_utils import save_yaml
from alphabase.constants.modification import MOD_DF

from peptdeep.settings import global_settings
from peptdeep.webui.ui_utils import (
    get_posix, select_files, file_type_selectbox
)

from peptdeep.webui.server import queue_folder

def nce_search():
    grid_nce_first = st.number_input(label='Start NCE for grid NCE search',value = global_settings['model_mgr']['transfer']['grid_nce_first']*1.0,step = 1.0)
    global_settings['model_mgr']['transfer']['grid_nce_first'] = grid_nce_first
    grid_nce_last = st.number_input(label='End NCE for grid NCE search',min_value = grid_nce_first, value = global_settings['model_mgr']['transfer']['grid_nce_last']*1.0,step = 1.0)
    global_settings['model_mgr']['transfer']['grid_nce_last'] = grid_nce_last
    grid_nce_step = st.number_input(label='Step NCE for grid NCE search', value = global_settings['model_mgr']['transfer']['grid_nce_step']*1.0,step = 1.0)
    global_settings['model_mgr']['transfer']['grid_nce_step'] = grid_nce_step

    grid_instrument = st.multiselect(label='Instruments for grid NCE search',
        options=global_settings['model_mgr']['instrument_group'],
        default = global_settings['model_mgr']['transfer']['grid_instrument']) 
    global_settings['model_mgr']['transfer']['grid_instrument'] = grid_instrument

def fine_tune():
    epoch_ms2 = st.number_input(label='Epoch to train MS2 model', value = global_settings['model_mgr']['transfer']['epoch_ms2'])
    global_settings['model_mgr']['transfer']['epoch_ms2'] = epoch_ms2
    warmup_epoch_ms2 = st.number_input(label='Warmup epoch to train MS2 model', value = global_settings['model_mgr']['transfer']['warmup_epoch_ms2'], max_value=epoch_ms2)
    global_settings['model_mgr']['transfer']['warmup_epoch_ms2'] = warmup_epoch_ms2
    batch_size_ms2 = st.number_input(label='Mini-batch size to train MS2 model', value = global_settings['model_mgr']['transfer']['batch_size_ms2'])
    global_settings['model_mgr']['transfer']['batch_size_ms2'] = batch_size_ms2
    lr_ms2 = st.number_input(label='Learning rate to train MS2 model', value = global_settings['model_mgr']['transfer']['lr_ms2'], format='%e', step=1e-5)
    global_settings['model_mgr']['transfer']['lr_ms2'] = lr_ms2
    
    epoch_rt_ccs = st.number_input(label='Epoch to train RT and CCS models', value = global_settings['model_mgr']['transfer']['epoch_rt_ccs'])
    global_settings['model_mgr']['transfer']['epoch_rt_ccs'] = epoch_rt_ccs
    warmup_epoch_rt_ccs = st.number_input(label='Warmup epoch to train RT and CCS model', value = global_settings['model_mgr']['transfer']['warmup_epoch_rt_ccs'], max_value=epoch_rt_ccs)
    global_settings['model_mgr']['transfer']['warmup_epoch_rt_ccs'] = warmup_epoch_rt_ccs
    batch_size_rt_ccs = st.number_input(label='Mini-batch size to train RT and CCS model', value = global_settings['model_mgr']['transfer']['batch_size_rt_ccs'])
    global_settings['model_mgr']['transfer']['batch_size_rt_ccs'] = batch_size_rt_ccs
    lr_rt_ccs = st.number_input(label='Learning rate to train RT and CCS model', value = global_settings['model_mgr']['transfer']['lr_rt_ccs'], format='%e', step=1e-5)
    global_settings['model_mgr']['transfer']['lr_rt_ccs'] = lr_rt_ccs

def add_other_psm_reader_mods():
    st.write("#### Other modification mapping for PSM readers")
    st.write('Thus PeptDeep supports to any modifications from other PSM readers')
    other_mod_expander = st.expander(label="Add other modification mapping")
    with other_mod_expander:
        mod_name = st.selectbox(label='AlphaBase modification',
            options=MOD_DF.index.values,
        )
        other_mods = st.text_input(
            label='Other modifications, sep by ";" for multiple ones',
            key='other_reader_mods'
        ).strip()
        st.text("Examples of other modifications: _(Dimethyl-n-0);_(Dimethyl) or K(Dimethyl-K-0)")

        if st.button("Add a modification mapping"):
            global_settings['model_mgr']['transfer'][
                'other_modification_mapping'
            ][mod_name] = other_mods.split(';')

        st.dataframe(pd.DataFrame().from_dict(
            global_settings['model_mgr']['transfer'][
                'other_modification_mapping'
            ],
            orient = 'index',
        ))

        def _clear_user_mods():
            global_settings['model_mgr']['transfer'][
                'other_modification_mapping'
            ] = {}
            st.session_state.other_reader_mods = ''

        st.button(label='Clear all other modification mapping', 
            on_click=_clear_user_mods
        )

def show():
    st.write("# Transfer learning")

    model_output_folder = st.text_input(
        label='Model output folder', 
        value=global_settings['model_mgr']['transfer']['model_output_folder'].format(
            PEPTDEEP_HOME=global_settings['PEPTDEEP_HOME']
        )
    )
    model_output_folder = get_posix(
        model_output_folder
    )
    global_settings['model_mgr']['transfer']['model_output_folder'] = model_output_folder

    st.write("### PSM files for training")
    psm_type = file_type_selectbox(
        ui_label='PSM type',
        st_key='select_psm_type',
        default_type=global_settings['model_mgr']['transfer']['psm_type'],
        monitor_files=global_settings['model_mgr']['transfer']['psm_files'],
        choices=global_settings['model_mgr']['transfer']['psm_type_choices'],
        index=global_settings['model_mgr']['transfer']['psm_type_choices'].index(
            global_settings['model_mgr']['transfer']['psm_type']
        )
    )
    global_settings['model_mgr']['transfer']['psm_type'] = psm_type

    psm_type_to_ext_dict = {
        "alphapept": ".ms_data.hdf",
        "pfind": ".spectra",
        "maxquant": "msms.txt",
        "diann": "tsv",
        "speclib_tsv": "tsv",
    }
    global_settings['model_mgr']['transfer']['psm_type'] = psm_type
    select_files(
        global_settings['model_mgr']['transfer']['psm_files'], 
        psm_type_to_ext_dict[psm_type], 
        "Input PSM files"
    )
    add_other_psm_reader_mods()

    st.write("### MS files for training")
    ms_file_type = file_type_selectbox(
        ui_label='MS file type',
        st_key='select_ms_file_type',
        default_type=global_settings['model_mgr']['transfer']['ms_file_type'],
        monitor_files=global_settings['model_mgr']['transfer']['ms_files'],
        choices=global_settings['model_mgr']['transfer']['ms_file_type_choices'], 
        index=global_settings['model_mgr']['transfer']['ms_file_type_choices'].index(
            global_settings['model_mgr']['transfer']['ms_file_type']
        )
    )
    ms_type_to_ext_dict = {
        "alphapept_hdf": ".ms_data.hdf",
        "thermo_raw": ".raw",
        "mgf": ".mgf",
        "mzml": ".mzml",
    }
    global_settings['model_mgr']['transfer']['ms_file_type'] = ms_file_type
    select_files(
        global_settings['model_mgr']['transfer']['ms_files'], 
        ms_type_to_ext_dict[ms_file_type], 
        "Input MS files"
    )

    st.write("### Training settings")

    training_expander = st.expander("Training hyper-parameters")
    with training_expander:
        fine_tune()

        global_settings['model_mgr']['transfer'][
            'psm_num_to_train_ms2'
        ] = st.number_input(
            label='PSM num to refine MS2 model', 
            value = int(global_settings['model_mgr']['transfer'][
                'psm_num_to_train_ms2'
            ]), step = 1
        )
        global_settings['model_mgr']['transfer'][
            'psm_num_per_mod_to_train_ms2'
        ] = st.number_input(
            label='PSM num per mod to refine MS2 model', 
            value = int(global_settings['model_mgr']['transfer'][
                'psm_num_per_mod_to_train_ms2'
            ]), step = 1
        )
        global_settings['model_mgr']['transfer'][
            'psm_num_to_test_ms2'
        ] = st.number_input(
            label='PSM num to test MS2 model', 
            value = int(global_settings['model_mgr']['transfer'][
                'psm_num_to_test_ms2'
            ]), step = 1
        )

        global_settings['model_mgr']['transfer'][
            'psm_num_to_train_rt_ccs'
        ] = st.number_input(
            label='PSM num to refine RT and CCS model', 
            value = int(global_settings['model_mgr']['transfer'][
                'psm_num_to_train_rt_ccs'
            ]), step = 1
        )
        global_settings['model_mgr']['transfer'][
            'psm_num_per_mod_to_train_rt_ccs'
        ] = st.number_input(
            label='PSM num per mod to refine RT and CCS model', 
            value = int(global_settings['model_mgr']['transfer'][
                'psm_num_per_mod_to_train_rt_ccs'
            ]), step = 1
        )
        
        global_settings['model_mgr']['transfer'][
            'psm_num_to_test_rt_ccs'
        ] = st.number_input(
            label='PSM num to test RT and CCS model', 
            value = int(global_settings['model_mgr']['transfer'][
                'psm_num_to_test_rt_ccs'
            ]), step = 1
        )
        global_settings['model_mgr']['transfer'][
            'top_n_mods_to_train'
        ] = st.number_input(
            label='Top n mods to refine models', 
            value = int(global_settings['model_mgr']['transfer'][
                'top_n_mods_to_train'
            ]), step = 1
        )

    st.write('#### Grid search for NCEs and instruments')
    st.write('If NCE and instrument are unknown, grid search will look for the best NCE and instrument)')
    grid_nce_search = bool(st.checkbox(label='Enabled', 
        value=global_settings['model_mgr']['transfer']['grid_nce_search']
    ))
    global_settings['model_mgr']['transfer']['grid_nce_search'] = grid_nce_search
    if grid_nce_search is True:
        nce_search()

    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d--%H-%M-%S.%f")
    task_name = st.text_input(label="Task name", value=f"peptdeep_transfer_{current_time}")
    
    if st.button(label='Submit for transfer learning'):
        global_settings['task_type'] = 'train'

        if not os.path.exists(model_output_folder):
            os.makedirs(model_output_folder)

        yaml_path = f'{queue_folder}/{task_name}.yaml'
        save_yaml(
            yaml_path, global_settings
        )
        save_yaml(
            os.path.join(
                model_output_folder, 
                f'{task_name}.yaml'
            ), 
            global_settings
        )
        
        st.write(f'`train` task saved as "{os.path.expanduser(yaml_path)}" in the task queue')