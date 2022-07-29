import streamlit as st
import pandas as pd
import os
from peptdeep.settings import global_settings
from peptdeep.cli import generate_library
from alphabase.constants.modification import MOD_DF

from peptdeep.webui.ui_utils import (
    files_in_pandas, update_input_paths
)

def mod_options():
    with st.form("Select modifications"):
        fixmod, = st.multiselect(
                'Please select fixed modifications',
                MOD_DF.index.values,
                default = global_settings['library']['input']['fix_mods']
            ),
        varmod, = st.multiselect(
                'Please select variable modifications',
                MOD_DF.index.values,
                default = global_settings['library']['input']['var_mods']
            ),
        global_settings['library']['input']['fix_mods'] = fixmod
        global_settings['library']['input']['var_mods'] = varmod

        st.form_submit_button("Click to confirm the selected modifications")
        st.text("Selected modifications")
        st.dataframe(MOD_DF.loc[fixmod+varmod,['mod_name','classification','composition','mass','modloss_composition','modloss','modloss_importance']])

        global_settings['common']['modloss_importance_level'] = st.number_input(
            'Modification loss importance level (for a PTM, fragment modloss mz=0 if modloss_importance<modloss_importance_level)', 
            value=global_settings['common']['modloss_importance_level']
        )

def varmod_range():
    min_varmod = st.number_input('Min number of variable modifications',min_value = 0, max_value = 1, value = global_settings['library']['input']['min_var_mod_num'], step = 1)
    max_varmod = st.number_input('Max number of variable modifications',min_value = 1, max_value = 2, value = global_settings['library']['input']['max_var_mod_num'], step = 1)
    global_settings['library']['input']['min_var_mod_num'] = min_varmod
    global_settings['library']['input']['max_var_mod_num'] = max_varmod


def choose_precursor_charge():
    from_charge = st.number_input('Min precursor charge', min_value = 1, max_value = 4, value = global_settings['library']['input']['min_precursor_charge'], step = 1)
    to_charge = st.number_input(
        'Max precursor charge',
        min_value = from_charge, max_value = 7, value = global_settings['library']['input']['max_precursor_charge'], step = 1
    )
    global_settings['library']['input']['min_precursor_charge'] = from_charge
    global_settings['library']['input']['max_precursor_charge'] = to_charge

def choose_precursor_mz():
    min_precursor_mz = st.number_input('Min precursor mz', value = global_settings['library']['input']['min_precursor_mz'])
    global_settings['library']['input']['min_precursor_mz'] = min_precursor_mz
    max_precursor_mz = st.number_input('Max precursor mz', min_value = min_precursor_mz, value = global_settings['library']['input']['max_precursor_mz'])
    global_settings['library']['input']['max_precursor_mz'] = max_precursor_mz

def add_decoy():
    decoy = st.selectbox('Decoy method',global_settings['library']['input']['decoy_choices'],index = 0)
    global_settings['library']['input']['decoy'] = decoy

def choose_protease():
    protease = st.selectbox(
        'Protease',
        global_settings['library']['input']['fasta']['protease_choices'],
    )
    global_settings['library']['input']['fasta']['protease'] = protease
    max_miss_cleave = st.number_input('Max number of miss cleavages',value = global_settings['library']['input']['fasta']['max_miss_cleave'])
    global_settings['library']['input']['fasta']['max_miss_cleave'] = max_miss_cleave

def choose_peptide_len():
    min_peptide_len = st.number_input('Min peptide length', value = global_settings['library']['input']['min_peptide_len'])
    max_peptide_len = st.number_input('Max peptide length', min_value = min_peptide_len, value = global_settings['library']['input']['max_peptide_len'])
    global_settings['library']['input']['min_peptide_len'] = min_peptide_len
    global_settings['library']['input']['max_peptide_len'] = max_peptide_len

def choose_frag_types():
    frag_types = st.multiselect(
        'Fragment types',(global_settings['model']['frag_types']),
        default = global_settings['library']['input']['frag_types']
    )
    global_settings['library']['input']['frag_types'] = frag_types
    max_frag_charge = st.number_input('Max fragment charge',min_value = 1, max_value = 2, value = global_settings['library']['input']['max_frag_charge'], step = 1)
    global_settings['library']['input']['max_frag_charge'] = max_frag_charge

def output_tsv():
    min_fragment_mz = st.number_input('Min fragment mz:', value = global_settings['library']['output_tsv']['min_fragment_mz'])
    global_settings['library']['output_tsv']['min_fragment_mz'] = min_fragment_mz
    max_fragment_mz = st.number_input('Max fragment mz:', min_value = min_fragment_mz, value = global_settings['library']['output_tsv']['max_fragment_mz'])
    global_settings['library']['output_tsv']['max_fragment_mz'] = max_fragment_mz
    min_relative_intensity = st.number_input('Min relative intensity:', value = global_settings['library']['output_tsv']['min_relative_intensity'])
    global_settings['library']['output_tsv']['min_relative_intensity'] = min_relative_intensity
    keep_higest_k_peaks = st.number_input('Number of highest peaks to keep:', value = global_settings['library']['output_tsv']['keep_higest_k_peaks'])
    global_settings['library']['output_tsv']['keep_higest_k_peaks'] = keep_higest_k_peaks
    global_settings['library']['output_tsv']['translate_mod_to_unimod_id']=bool(st.checkbox('Translate modifications to Unimod ids'))

def select_files(_input_type):
    path = st.text_input(f"File paths ({_input_type if _input_type=='fasta' else 'tsv/csv/txt'} files)")
    col1, col2, col3 = st.columns([0.5,0.5,2])
    with col1:
        add = st.button('Add')
    with col2:
        remove = st.button('Remove')
    with col3:
        clear = st.button('Clear all files')
    if add is True:
        if path not in global_settings['library']['input']['paths']:
            global_settings['library']['input']['paths'].append(path)
    if remove is True:
        if path in global_settings['library']['input']['paths']:
            global_settings['library']['input']['paths'].remove(path)
    if clear is True:
        global_settings['library']['input']['paths'] = []
    update_input_paths(global_settings['library']['input']['paths'])
    st.table(files_in_pandas(global_settings['library']['input']['paths']))

def input_files():
    def on_input_change():
        if len(global_settings['library']['input']['paths'])>0:
            st.warning("Please clear all input files before changing the input file type")
            st.session_state.input_type = global_settings['library']['input']['type']
            
    update_input_paths(global_settings['library']['input']['paths'])
    input_type = st.selectbox(
        'Input file type',
        global_settings['library']['input']['type_choices'],
        key='input_type',
        on_change=on_input_change
    )
    global_settings['library']['input']['type'] = input_type

    select_files(input_type)
    return input_type

def show():
    st.write("# Library Prediction")

    st.write('### Input')

    _input_type = input_files()

    add_decoy()

    choose_precursor_mz()

    if _input_type == 'fasta':
        choose_protease()
        mod_options()
        varmod_range()
        choose_precursor_charge()
        choose_peptide_len()

    elif _input_type == 'sequence_table':
        mod_options()
        varmod_range()
        choose_precursor_charge()

    elif _input_type == 'peptide_table':
        choose_precursor_charge()

    choose_frag_types()


    st.write("### Output")

    output_folder = st.text_input("Output folder", value=global_settings['library']['output_folder'])
    output_folder = os.path.abspath(os.path.expanduser(os.path.expandvars(output_folder)))
    global_settings['library']['output_folder'] = output_folder

    tsv_enabled = bool(st.checkbox('Output TSV (for DiaNN/Spectronaut)'))
    global_settings['library']['output_tsv']['enabled'] = tsv_enabled
    if tsv_enabled:
        output_tsv()

    if st.button('Generate library'):
        if len(global_settings['library']['input']['paths']) > 0:
            generate_library()
            st.write('Library generated!')
        else:
            st.warning(f'Please select the input {_input_type} files')
