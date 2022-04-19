from email.policy import default
from multiprocessing.sharedctypes import Value
import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep
import pandas as pd
import os
import time
from peptdeep.settings import global_settings

def mod_options():
    fixmod = st.multiselect(
            'Please select fix_mod',
            ['Carbamidomethyl@C', 'B', 'C', 'D'],
            default = ['Carbamidomethyl@C']
        ), 
    varmod = st.multiselect(
            'Please select var_mod',
            ['Oxidation@M', 'B', 'C', 'D'],
            default = ['Oxidation@M']
        ), 
    global_settings['library']['input']['fix_mods'] = fixmod
    global_settings['library']['input']['var_mods'] = varmod


def varmod_range():
    min_varmod = st.number_input('min_var_mod_num',min_value = 0, max_value = 1, value = 0, step = 1)
    max_varmod = st.number_input('max_var_mod_num',min_value = 1, max_value = 2, value = 2, step = 1)
    global_settings['library']['input']['min_var_mod_num'] = min_varmod
    global_settings['library']['input']['max_var_mod_num'] = max_varmod


def choose_precursor_charge():
    from_charge = st.number_input('precursor charge from',min_value = 1, max_value = 4, value = 2, step = 1)
    to_charge = st.number_input(
        'to',
        min_value = from_charge, max_value = 6, value = from_charge, step = 1 
    )
    global_settings['library']['input']['min_precursor_charge'] = from_charge
    global_settings['library']['input']['max_precursor_charge'] = to_charge

def choose_precursor_mz():
    min_precursor_mz = st.number_input('min precursor mz', value = 400)
    global_settings['library']['input']['min_precursor_mz'] = min_precursor_mz
    max_precursor_mz = st.number_input('max precursor mz', min_value = min_precursor_mz, value = 2000)
    global_settings['library']['input']['max_precursor_mz'] = max_precursor_mz

def add_decoy():
    add_decoy = st.checkbox('Add decoy')
    if add_decoy == 1:
        decoy = st.selectbox('decoy',global_settings['library']['input']['decoy_choices'],index = 2)
        global_settings['library']['input']['decoy'] = decoy

def choose_protease():
    protease = st.selectbox(
        'protease', 
        global_settings['library']['input']['fasta']['protease_choices'],
    )
    global_settings['library']['input']['fasta']['protease'] = protease
    max_miss_cleave = st.number_input('max_miss_cleave',value = 2)
    global_settings['library']['input']['fasta']['max_miss_cleave'] = max_miss_cleave

def choose_peptide_len():
    min_peptide_len = st.number_input('min peptide len:', value = 7)
    max_peptide_len = st.number_input('max peptide len:', min_value = min_peptide_len, value = 35)
    global_settings['library']['input']['min_peptide_len'] = min_peptide_len
    global_settings['library']['input']['max_peptide_len'] = max_peptide_len

def choose_frag_types():
    frag_types = st.multiselect(
        'frag_types',(global_settings['model']['frag_types']),
        default = ['b','y']
    )
    global_settings['library']['input']['frag_types'] = frag_types
    max_frag_charge = st.number_input('max fragment charge:',min_value = 1, max_value = 2, value = 2, step = 1)
    global_settings['library']['input']['max_frag_charge'] = max_frag_charge

def output_tsv():
    min_fragment_mz = st.number_input('min fragment mz:', value = 300)
    global_settings['library']['output_tsv']['min_fragment_mz'] = min_fragment_mz
    max_fragment_mz = st.number_input('max fragment mz:', min_value = min_fragment_mz, value = 2000)
    global_settings['library']['output_tsv']['max_fragment_mz'] = max_fragment_mz
    min_relative_intensity = st.number_input('min relative intensity:', value = 0.02)
    global_settings['library']['output_tsv']['min_relative_intensity'] = min_relative_intensity
    keep_higest_k_peaks = st.number_input('top n peaks:', value = 12)
    global_settings['library']['output_tsv']['keep_higest_k_peaks'] = keep_higest_k_peaks

def input_type():
    input_type = st.selectbox(
        'input types',
        global_settings['library']['input']['type_choices'],
    )
    global_settings['library']['input']['type'] = input_type
    return input_type

def show():

    st.write("# Library started")
    st.text("Welcome to library of AlphaPeptDeep.")

    st.write('### input')

    _input_type = input_type()

    path = st.text_input('paths')

    add_decoy()

    choose_precursor_mz()

    if _input_type == 'fasta':
        choose_protease()
        mod_options()
        varmod_range()
        choose_precursor_charge()
        choose_peptide_len()


    if _input_type == 'sequence_table':
        mod_options()
        varmod_range()
        choose_precursor_charge()

    if _input_type == 'peptide_table':
        choose_precursor_charge()
    
    choose_frag_types()


    st.write("### output")

    if st.checkbox('Output TSV') == 1:
        output_tsv()
    
