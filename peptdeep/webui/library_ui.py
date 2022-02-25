import streamlit as st
from peptdeep.webui.ui_utils import markdown_link
import peptdeep
import pandas as pd
import os
import time
from peptdeep.settings import global_settings

def show():

    st.write("# Library started")
    st.text("Welcome to library of AlphaPeptDeep.")
    st.write('### input')

    input_type = st.selectbox(
        'input types',
        global_settings['library']['input']['type_choices'],
    )
    #st.write('You selected:', input_type)
    global_settings['library']['input']['type'] = input_type

    path = st.text_input('paths')
    add_decoy = st.checkbox('Add decoy')
    decoy_choices = st.selectbox('decoy_choices',['pseudo_reverse','diann'])

    if input_type == 'fasta':
        enzyme = st.text_input('enzyme')
        max_miss_cleave = st.text_input('max_miss_cleave')
        fixmod_options = st.multiselect(
            'Please select fixmod',
            ['Carbamidomethyl@C', 'B', 'C', 'D'],
            [])
        global_settings['library']['input']['fixmod'] = fixmod_options
        #st.write('You selected:', fixmod_options)
        varmod_options = st.multiselect(
            'Please select varmod',
            ['Oxidation@M','test1','test2','test3'],
            []
        )
        global_settings['library']['input']['varmod'] = varmod_options

        min_varmod = st.selectbox(
            'min_varmod',
            ('0','1')
        )
        global_settings['library']['input']['min_varmod'] = min_varmod
        max_varmod = st.selectbox(
            'max_varmod',
            ('1','2')
        )
        global_settings['library']['input']['max_varmod'] = max_varmod


        from_charge = st.selectbox('precursor charge from',
        ('1','2','3','4'))
        to_charge = st.selectbox(
            'to',
            ('1','2','3','4','5','6')
        )
        if to_charge < from_charge:
            st.text("ERROR, num of to_charge should be larger than from_charge")
            st.text('Please select again.')
        global_settings['library']['output']['min_precursor_charge'] = from_charge
        global_settings['library']['output']['max_precursor_charge'] = to_charge
    if input_type == 'sequence_list':
        fixmod_options = st.multiselect(
            'Please select fixmod',
            ['Carbamidomethyl@C', 'B', 'C', 'D'],
            [])
        global_settings['library']['input']['fixmod'] = fixmod_options
    #st.write('You selected:', fixmod_options)

        varmod_options = st.multiselect(
            'Please select varmod',
            ['Oxidation@M','test1','test2','test3'],
            []
        )
        global_settings['library']['input']['varmod'] = varmod_options

        min_varmod = st.selectbox(
            'min_varmod',
            ('0','1')
        )
        global_settings['library']['input']['min_varmod'] = min_varmod
        max_varmod = st.selectbox(
            'max_varmod',
            ('1','2')
        )
        global_settings['library']['input']['max_varmod'] = max_varmod


        from_charge = st.selectbox('precursor charge from',
        ('1','2','3','4'))
        to_charge = st.selectbox(
            'to',
            ('1','2','3','4','5','6')
        )
        if to_charge < from_charge:
            st.text("ERROR, num of to_charge should be larger than from_charge")
            st.text('Please select again.')
        global_settings['library']['output']['min_precursor_charge'] = from_charge
        global_settings['library']['output']['max_precursor_charge'] = to_charge
    if input_type == 'peptide_list':
        from_charge = st.selectbox('precursor charge from',
        ('1','2','3','4'))
        to_charge = st.selectbox(
            'to',
            ('1','2','3','4','5','6')
        )
        if to_charge < from_charge:
            st.text("ERROR, num of to_charge should be larger than from_charge")
            st.text('Please select again.')
        global_settings['library']['output']['min_precursor_charge'] = from_charge
        global_settings['library']['output']['max_precursor_charge'] = to_charge
   




    st.write("### output")

    output_type = st.selectbox(
        'output type',
        ('hdf','TSV(spectronant/diann)')
    )
    global_settings['library']['output']['type'] = output_type
