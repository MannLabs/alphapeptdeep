# from alphapept.gui

import streamlit as st

from PIL import Image
import os
import socket
import peptdeep

from peptdeep.webui import (
    model_ui, startpage, rescore_ui, 
    library_ui, transfer_ui,
    settings_ui, server_ui,
)

_this_file = __file__
_this_directory = os.path.dirname(_this_file)
LOGO_PATH = os.path.join(_this_directory, 'logos', 'peptdeep.png')
ICON_PATH = os.path.join(_this_directory, 'logos', 'peptdeep.ico')
image = Image.open(LOGO_PATH)
icon = Image.open(ICON_PATH)
computer_name = socket.gethostname()

st.set_page_config(
    page_title=f"PeptDeep {peptdeep.__version__}",
    page_icon=icon,
    layout="wide",
)

hide_streamlit_style = """
<style>
#MainMenu {visibility: hidden;}
footer {visibility: hidden;}
</style>

"""
st.markdown(hide_streamlit_style, unsafe_allow_html=True)

st.sidebar.image(image, width = 300)
st.sidebar.code(f"AlphaPeptDeep (PeptDeep) {peptdeep.__version__} \n{computer_name}")

sidebar = {
    'Start page': startpage.show,
    'Model': model_ui.show,
    'Transfer': transfer_ui.show,
    'Library': library_ui.show,
    'Rescore': rescore_ui.show,
    'Server': server_ui.show,
    'Settings':  settings_ui.show,
}

menu = st.sidebar.radio("", list(sidebar.keys()))

if menu:
    sidebar[menu]()

link = f'[PeptDeep on GitHub]({peptdeep.__github__})'
st.sidebar.markdown(link, unsafe_allow_html=True)