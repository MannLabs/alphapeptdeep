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
    min_varmod = st.selectbox(
            'min_var_mod_num',
            ('0','1'),
            index = 0
        ),
    max_varmod = st.selectbox(
            'max_var_mod_num',
            ('1','2'),
            index = 1
        ),
    global_settings['library']['input']['min_var_mod_num'] = min_varmod
    global_settings['library']['input']['max_var_mod_num'] = max_varmod


def choose_precursor_charge():
    from_charge = st.selectbox('precursor charge from',
    ('1','2','3','4'),index = 1)
    options = range(int(from_charge),int(7))
    to_charge = st.selectbox(
        'to',
        (options),index = 2
    )
    global_settings['library']['output']['min_precursor_charge'] = from_charge
    global_settings['library']['output']['max_precursor_charge'] = to_charge

def choose_precursor_mz():
    min_precursor_mz = st.text_input('min precursor mz', value = "400")
    global_settings['library']['input']['min_precursor_mz'] = min_precursor_mz
    max_precursor_mz = st.text_input('max precursor mz', value = "2000")
    global_settings['library']['input']['max_precursor_mz'] = max_precursor_mz

def add_decoy():
    add_decoy = st.checkbox('Add decoy')
    if add_decoy == 1:
        decoy = st.selectbox('decoy',
        global_settings['library']['input']['decoy_choices']
        )
        global_settings['library']['input']['decoy'] = decoy

def choose_protease():
    protease = st.text_input('protease', value = 'trypsin')
    global_settings['library']['input']['fasta']['protease'] = protease
    max_miss_cleave = st.text_input('max_miss_cleave',value = '2')
    global_settings['library']['input']['fasta']['max_miss_cleave'] = max_miss_cleave

def choose_peptide_len():
    min_peptide_len = st.text_input('min peptide len:', value = 7)
    max_peptide_len = st.text_input('max peptide len:', value = 35)
    global_settings['library']['input']['min_peptide_len'] = min_peptide_len
    global_settings['library']['input']['max_peptide_len'] = max_peptide_len

def choose_frag_types():
    frag_types = st.multiselect(
        'frag_types',('b','y','b_modloss','y_modloss'),
        default = ['b','y']
    )
    global_settings['library']['input']['frag_types'] = frag_types
    max_frag_charge = st.selectbox('max fragment charge:',['1','2'], index = 1)
    global_settings['library']['input']['max_frag_charge'] = max_frag_charge

def output_tsv():
    min_fragment_mz = st.text_input('min fragment mz:', value = "300")
    global_settings['library']['output']['min_fragment_mz'] = min_fragment_mz
    max_fragment_mz = st.text_input('max fragment mz:', value = "2000")
    global_settings['library']['output']['max_fragment_mz'] = max_fragment_mz
    min_relative_intensity = st.text_input('min relative intensity:', value = "0.02")
    global_settings['library']['output']['min_relative_intensity'] = min_relative_intensity
    keep_higest_k_peaks = st.text_input('top n peaks:', value = "12")
    global_settings['library']['output']['keep_higest_k_peaks'] = keep_higest_k_peaks

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
    
