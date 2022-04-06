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
    return (
        st.multiselect(
            'Please select fixmod',
            ['Carbamidomethyl@C', 'B', 'C', 'D'],
            default = ['Carbamidomethyl@C']
        ), 
        st.multiselect(
            'Please select varmod',
            ['Oxidation@M', 'B', 'C', 'D'],
            default = ['Oxidation@M']
        ), 
    )

def varmod_range():
    return(
        st.selectbox(
            'min_varmod',
            ('0','1'),
            index = 0
        ),
        st.selectbox(
            'max_varmod',
            ('1','2'),
            index = 1
        ),
    )



def show():

    st.write("# Library started")
    st.text("Welcome to library of AlphaPeptDeep.")
    st.write('### input')

    input_type = st.selectbox(
        'input types',
        global_settings['library']['input']['type_choices'],
    )
    global_settings['library']['input']['type'] = input_type

    path = st.text_input('paths')
  
    add_decoy = st.checkbox('Add decoy')
    if add_decoy == 1:
        decoy = st.selectbox('decoy',['pseudo_reverse','diann'])
        global_settings['library']['input']['decoy'] = decoy

    if input_type == 'fasta':
        enzyme = st.text_input('enzyme', value = 'trypsin')
        global_settings['library']['input']['fasta']['enzyme'] = enzyme
        max_miss_cleave = st.text_input('max_miss_cleave',value = '2')
        global_settings['library']['input']['fasta']['max_miss_cleave'] = max_miss_cleave
        fixmod, varmod = mod_options()
        global_settings['library']['input']['fixmod'] = fixmod
        global_settings['library']['input']['varmod'] = varmod

        min_varmod,max_varmod = varmod_range()
        global_settings['library']['input']['min_varmod'] = min_varmod
        global_settings['library']['input']['max_varmod'] = max_varmod


        from_charge = st.selectbox('precursor charge from',
        ('1','2','3','4'),index = 1)
        options = range(int(from_charge),int(7))
        to_charge = st.selectbox(
            'to',
            (options)
        )
        global_settings['library']['output']['min_precursor_charge'] = from_charge
        global_settings['library']['output']['max_precursor_charge'] = to_charge

        min_peptide_len = st.text_input('min peptide len:', value = 7)
        max_peptide_len = st.text_input('max peptide len:', value = 35)
        global_settings['library']['input']['min_peptide_len'] = min_peptide_len
        global_settings['library']['input']['max_peptide_len'] = max_peptide_len


    if input_type == 'sequence_table':
        fixmod, varmod = mod_options()
        global_settings['library']['input']['fixmod'] = fixmod
        global_settings['library']['input']['varmod'] = varmod

        min_varmod,max_varmod = varmod_range()
        global_settings['library']['input']['min_varmod'] = min_varmod
        global_settings['library']['input']['max_varmod'] = max_varmod

        from_charge = st.selectbox('precursor charge from',
        ('1','2','3','4'),index = 1)
        options = range(int(from_charge),int(7))
        to_charge = st.selectbox(
            'to',
            (options)
        )
        global_settings['library']['output']['min_precursor_charge'] = from_charge
        global_settings['library']['output']['max_precursor_charge'] = to_charge

    if input_type == 'peptide_table':
   
        from_charge = st.selectbox('precursor charge from',
        ('1','2','3','4'),index = 1)
        options = range(int(from_charge),int(7))
        to_charge = st.selectbox(
            'to',
            (options)
        )

        global_settings['library']['output']['min_precursor_charge'] = from_charge
        global_settings['library']['output']['max_precursor_charge'] = to_charge
   
    frag_types = st.multiselect(
        'frag_types',('b','y','b_modloss','y_modloss'),
        default = ['b','y']
    )
    global_settings['library']['input']['frag_types'] = frag_types
    max_frag_charge = st.selectbox('max fragment charge:',['1','2'], index = 1)
    global_settings['library']['input']['max_frag_charge'] = max_frag_charge


    st.write("### output")

    output_type = st.selectbox(
        'output type',
        ('hdf','TSV (spectronant/diann)')
    )
    global_settings['library']['output']['type'] = output_type

    if output_type == 'TSV (spectronant/diann)':
        min_fragment_mz = st.text_input('min fragment mz:', value = "300")
        global_settings['library']['output']['min_fragment_mz'] = min_fragment_mz
        max_fragment_mz = st.text_input('max fragment mz:', value = "2000")
        global_settings['library']['output']['max_fragment_mz'] = max_fragment_mz
        min_relative_intensity = st.text_input('min relative intensity:', value = "0.02")
        global_settings['library']['output']['min_relative_intensity'] = min_relative_intensity
        top_n_peaks = st.text_input('top n peaks:', value = "12")
        global_settings['library']['output']['top_n_peaks'] = top_n_peaks
        min_precursor_mz = st.text_input('min precursor mz', value = "400")
        global_settings['library']['output']['min_precursor_mz'] = min_precursor_mz
        max_precursor_mz = st.text_input('max precursor mz', value = "2000")
        global_settings['library']['output']['max_precursor_mz'] = max_precursor_mz

    
