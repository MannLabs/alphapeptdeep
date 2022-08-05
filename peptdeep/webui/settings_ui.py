import yaml
import streamlit as st
from io import StringIO

from peptdeep.settings import (
    global_settings,
    update_settings
)

def show():
    st.write("# Import previous yaml settings")
    uploaded_file = st.file_uploader("Choose a local yaml file")
    if uploaded_file is not None:
        f = StringIO(uploaded_file.getvalue().decode("utf-8"))
        uploaded_settings = yaml.load(f, Loader=yaml.FullLoader)
        update_settings(global_settings, uploaded_settings)
        st.write("Global settings have been updated")