import streamlit as st
import os
from alphabase.yaml_utils import load_yaml

from peptdeep.webui.server import get_yamls, queue_folder
from peptdeep.webui.ui_utils import files_in_pandas

from peptdeep.webui.server import _server

def display_tasks():
    st.write("## Tasks in the queue")

    yamls = get_yamls(queue_folder)
    st.write(f"There are {len(yamls)} tasks in the queue:")
    df = files_in_pandas(yamls)

    st.write(f"See `{os.path.expanduser(queue_folder)}`")

    tasks = []
    for _yml in yamls:
        _dict = load_yaml(_yml)
        if 'task_type' in _dict:
            task_type = _dict['task_type']
        else:
            task_type = 'library'
        tasks.append(task_type)

    df['Task Type'] = tasks

    st.table(df)

def show():
    st.write("# AlphaPeptDeep Server")

    if _server.process is not None:
        st.write("The server is running or waiting for tasks")
    else:
        st.write("The serve is not running")

    if st.button("Start the server"):
        _server.start()

    if st.button("Terminate the server"):
        _server.terminate()

    display_tasks()
