import streamlit as st
import pandas as pd
import os

from datetime import datetime

from alphabase.constants.modification import MOD_DF
from alphabase.yaml_utils import save_yaml

from peptdeep.settings import global_settings

from peptdeep.webui.ui_utils import (
    get_posix, select_files, file_type_selectbox
)

from peptdeep.webui.server import queue_folder

def mod_options():
    with st.form(key="Select modifications"):
        st.write('#### Fixed and variable modificatins')
        fixmod = st.multiselect(
            label='Please select fixed modifications',
            options=MOD_DF.index.values,
            default = global_settings['library']['input']['fix_mods']
        )
        varmod = st.multiselect(
            label='Please select variable modifications',
            options=MOD_DF.index.values,
            default = global_settings['library']['input']['var_mods']
        )
        global_settings['library']['input']['fix_mods'] = fixmod
        global_settings['library']['input']['var_mods'] = varmod

        st.form_submit_button(label="Click to add these selected modifications")
        st.write("Selected modifications:")
        st.dataframe(MOD_DF.loc[fixmod+varmod,[
            'mod_name','classification','composition','mass',
            'modloss_composition','modloss','modloss_importance'
        ]])

    varmod_range()

def varmod_range():
    min_varmod = st.number_input(label='Min number of variable modifications',
        value = global_settings['library']['input']['min_var_mod_num'], 
        min_value = 0, step = 1,
    )
    max_varmod = st.number_input(label='Max number of variable modifications',
        value = global_settings['library']['input']['max_var_mod_num'], 
        min_value = 0, step = 1,
    )
    global_settings['library']['input']['min_var_mod_num'] = min_varmod
    global_settings['library']['input']['max_var_mod_num'] = max_varmod

def specialmod_options():
    st.write('#### Special modificatins')
    st.write('*Useful for Phospho@S/T or GlyGly@K*')
    st.write('- For Phospho@S/T or HexNAc@S, as a sequence may generate many peptidoforms, this can control the overall number.')
    st.write('- For GlyGly@K or GG@K, it will not occur at C-term Lys/K, using `special modifications` to enable this feature.')
    specialmod_expander = st.expander(label='Special modificatins')
    with specialmod_expander:
        with st.form(key="Select special modifications"):
            global_settings['library']['input']['special_mods'] = st.multiselect(
                label='Please select special modifications',
                options=MOD_DF.index.values,
                default=global_settings['library']['input']['special_mods']
            )
            st.form_submit_button(label="Click to add selected modifications")
            st.write("Selected special modifications:")
            st.dataframe(MOD_DF.loc[
                global_settings['library']['input']['special_mods'],
                [
                    'mod_name','classification','composition','mass',
                    'modloss_composition','modloss','modloss_importance'
                ]
            ])

        specialmod_range()

def specialmod_range():
    min_specialmod = st.number_input(label='Min number of special modifications',
        value = global_settings['library']['input']['min_special_mod_num'], 
        min_value = 0, step = 1
    )
    max_specialmod = st.number_input(label='Max number of special modifications',
        value = global_settings['library']['input']['max_special_mod_num'], 
        min_value = 0, step = 1
    )
    global_settings['library']['input']['min_special_mod_num'] = min_specialmod
    global_settings['library']['input']['max_special_mod_num'] = max_specialmod

    st.write("Special modifications cannot modify AAs at:")
    st.write("*e.g. GlyGly@K will not occur at C-term Lys/K*")
    global_settings['library']['input'][
        'special_mods_cannot_modify_pep_n_term'
    ] = bool(
        st.checkbox(label='N-term', 
        value=global_settings['library']['input'][
            'special_mods_cannot_modify_pep_n_term'
        ])
    )
    global_settings['library']['input'][
        'special_mods_cannot_modify_pep_c_term'
    ] = bool(
        st.checkbox(label='C-term', 
        value=global_settings['library']['input'][
            'special_mods_cannot_modify_pep_c_term'
        ])
    )

def labeling_options():
    def _concat_df_dict(d):
        df_list = []
        for channel, mods in d.items():
            _df = MOD_DF.loc[mods,['mod_name','composition','mass']]
            _df['labeling_channel'] = channel
            df_list.append(_df)
        if len(df_list) == 0:
            return pd.DataFrame()
        else:
            return pd.concat(df_list, ignore_index=True)
    def _clear_all():
        global_settings['library']['input']['labeling_channels'] = {}
        st.session_state.select_labeling = []
        st.session_state.labeling_channel_id = ''
        return
    st.write('#### Peptide labeling')
    st.write('*For multiplex-DIA (mDIA) workflow*')
    labeling_expander = st.expander(label='Labeling channels')
    with labeling_expander:
        with st.form(key="Peptide labeling"):
            channel = st.text_input(label="Channel",key='labeling_channel_id')

            mods = st.multiselect(
                label='Please select labeling modifications',
                options=MOD_DF.index.values,
                key='select_labeling'
            )
            if channel and len(mods) > 0:
                try:
                    channel = int(channel)
                except ValueError:
                    pass
                global_settings['library']['input']['labeling_channels'][channel] = mods

            st.form_submit_button(label="Add selected labeling")
            st.write("Selected labeling modifications:")
            st.dataframe(_concat_df_dict(global_settings['library']['input']['labeling_channels']))
        st.button(label='Clear all labeling', on_click=_clear_all)

def choose_precursor_charge():
    from_charge = st.number_input(label='Min precursor charge', min_value = 1, max_value = 4, value = global_settings['library']['input']['min_precursor_charge'], step = 1)
    to_charge = st.number_input(
        label='Max precursor charge',
        min_value = from_charge, max_value = 7, value = global_settings['library']['input']['max_precursor_charge'], step = 1
    )
    global_settings['library']['input']['min_precursor_charge'] = from_charge
    global_settings['library']['input']['max_precursor_charge'] = to_charge

def choose_precursor_mz():
    min_precursor_mz = st.number_input(label='Min precursor mz', value = global_settings['library']['input']['min_precursor_mz'])
    global_settings['library']['input']['min_precursor_mz'] = min_precursor_mz
    max_precursor_mz = st.number_input(label='Max precursor mz', min_value = min_precursor_mz, value = global_settings['library']['input']['max_precursor_mz'])
    global_settings['library']['input']['max_precursor_mz'] = max_precursor_mz

def add_decoy():
    decoy = st.selectbox(label='Decoy method',options=global_settings['library']['input']['decoy_choices'],index = 0)
    global_settings['library']['input']['decoy'] = decoy

def choose_protease():
    protease = st.selectbox(
        label='Protease',
        options=global_settings['library']['input']['fasta']['protease_choices'],
    )
    global_settings['library']['input']['fasta']['protease'] = protease
    max_miss_cleave = st.number_input(label='Max number of miss cleavages',value = global_settings['library']['input']['fasta']['max_miss_cleave'])
    global_settings['library']['input']['fasta']['max_miss_cleave'] = max_miss_cleave

def choose_peptide_len():
    min_peptide_len = st.number_input(label='Min peptide length', value = global_settings['library']['input']['min_peptide_len'])
    max_peptide_len = st.number_input(label='Max peptide length', min_value = min_peptide_len, value = global_settings['library']['input']['max_peptide_len'])
    global_settings['library']['input']['min_peptide_len'] = min_peptide_len
    global_settings['library']['input']['max_peptide_len'] = max_peptide_len

def choose_frag_types():
    frag_types = st.multiselect(
        label='Fragment types',options=(global_settings['model']['frag_types']),
        default = global_settings['library']['input']['frag_types']
    )
    global_settings['library']['input']['frag_types'] = frag_types
    max_frag_charge = st.number_input(label='Max fragment charge',min_value = 1, max_value = 2, value = global_settings['library']['input']['max_frag_charge'], step = 1)
    global_settings['library']['input']['max_frag_charge'] = max_frag_charge

def output_tsv():
    min_fragment_mz = st.number_input(label='Min fragment mz:', value = global_settings['library']['output_tsv']['min_fragment_mz'])
    global_settings['library']['output_tsv']['min_fragment_mz'] = min_fragment_mz
    max_fragment_mz = st.number_input(label='Max fragment mz:', min_value = min_fragment_mz, value = global_settings['library']['output_tsv']['max_fragment_mz'])
    global_settings['library']['output_tsv']['max_fragment_mz'] = max_fragment_mz
    min_relative_intensity = st.number_input(label='Min relative intensity:', value = global_settings['library']['output_tsv']['min_relative_intensity'])
    global_settings['library']['output_tsv']['min_relative_intensity'] = min_relative_intensity
    keep_higest_k_peaks = st.number_input(label='Number of highest peaks to keep:', value = global_settings['library']['output_tsv']['keep_higest_k_peaks'])
    global_settings['library']['output_tsv']['keep_higest_k_peaks'] = keep_higest_k_peaks
    global_settings['library']['output_tsv']['translate_mod_to_unimod_id']=bool(st.checkbox(label='Translate modifications to Unimod ids'))

def show():
    st.write("# Library Prediction")

    st.write('### Input')

    infile_type = file_type_selectbox(
        ui_label='Input file type',
        st_key='lib_input_type',
        default_type=global_settings['library']['input']['infile_type'],
        monitor_files=global_settings['library']['input']['infiles'],
        choices=global_settings['library']['input']['infile_type_choices'], 
        index=global_settings['library']['input']['infile_type_choices'].index(
            global_settings['library']['input']['infile_type']
        )
    )
    global_settings['library']['input']['infile_type'] = infile_type

    if infile_type != 'fasta':
        import pandas as pd
        df = pd.DataFrame({
            'sequence': ['ACDEFGHIK','LMNPQRSTVK','WYVSTR'],
            'mods': ['Carbamidomethyl@C','Acetyl@Protein N-term;Phospho@S',''],
            'mod_sites': ['2','0;7',''],
            'charge': [2,3,1],
        })
        infile_expander = st.expander("Input file examples")
        with infile_expander:
            st.write('`sequence_table`:')
            st.dataframe(df[['sequence']])
            st.write('`peptide_table`:')
            st.dataframe(df[['sequence','mods','mod_sites']])
            st.write('`precursor_table`:')
            st.dataframe(df[['sequence','mods','mod_sites','charge']])

    infile_ext_dict = {
        'fasta': ['.fasta','.fa'],
        'sequence_table': ['tsv','txt','csv'],
        'peptide_table': ['tsv','txt','csv'],
        'precursor_table': ['tsv','txt','csv'],
    }
    select_files(
        global_settings['library']['input']['infiles'],
        infile_ext_dict[infile_type],
        'Input sequence files',
    )

    st.write('### Library settings')
    add_decoy()

    if infile_type == 'fasta':
        choose_protease()
        mod_options()
        specialmod_options()

    elif infile_type == 'sequence_table':
        mod_options()
        specialmod_options()
    
    labeling_options()

    st.write("#### Common peptide settings")
    
    if infile_type == 'fasta':
        choose_peptide_len()

    if infile_type in ['fasta','sequence_table','peptide_table']:
        choose_precursor_charge()
    
    choose_precursor_mz()
    choose_frag_types()

    st.write("### Output")

    output_folder = st.text_input(
        label="Output folder", 
        value=global_settings['library']['output_folder'].format(
            PEPTDEEP_HOME=global_settings['PEPTDEEP_HOME']
        )
    )
    output_folder = os.path.expanduser(output_folder)
    output_folder = get_posix(output_folder)
    global_settings['library']['output_folder'] = output_folder

    tsv_enabled = bool(st.checkbox(label='Output TSV (for DiaNN/Spectronaut)', value=global_settings['library']['output_tsv']['enabled']))
    st.warning("Writing the TSV file for a big library is very slow")
    global_settings['library']['output_tsv']['enabled'] = tsv_enabled
    if tsv_enabled:
        output_tsv()

    now = datetime.now()
    current_time = now.strftime("%Y-%m-%d--%H-%M-%S.%f")
    task_name = st.text_input(label="Task name", value=f"peptdeep_library_{current_time}")

    if st.button(label='Submit for library prediction'):
        global_settings['task_type'] = 'library'

        if not os.path.exists(output_folder):
            os.makedirs(output_folder)

        yaml_path = f'{queue_folder}/{task_name}.yaml'
        save_yaml(
            yaml_path, global_settings
        )
        save_yaml(
            os.path.join(
                output_folder, 
                f'{task_name}.yaml'
            ), 
            global_settings
        )
        st.write(f'`library` task saved as `{os.path.expanduser(yaml_path)}`')
