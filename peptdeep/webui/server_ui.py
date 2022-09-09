import streamlit as st
import os
import psutil
import time

from alphabase.yaml_utils import load_yaml

from peptdeep.webui.server import (
    get_yamls, queue_folder, _server
)
from peptdeep.webui.ui_utils import files_in_pandas

def display_tasks():
    st.write("## Tasks in the queue")

    yamls = get_yamls(queue_folder)
    st.write(f"There are {len(yamls)} tasks in the queue:")
    df = files_in_pandas(yamls)

    st.write(f"See yaml files in `{os.path.expanduser(queue_folder)}`")

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

    if st.button("Start the server"):
        _server.start()
        print("Server started")

    if st.button("Stop the server"):
        _server.terminate()
        print("Server stopped")

    if _server.process is not None:
        st.write("The server is running or waiting for tasks")
        st.warning("Stop the server before exit")
    else:
        st.write("The server is not running")
        st.warning("Start the server to submit tasks")

    display_tasks()

    st.write("### Hardware utilization")

    c1,c2 = st.columns(2)
    c1.text("Ram")
    ram = c1.progress(0)
    c2.text("CPU")
    cpu = c2.progress(0)

    while True:
        ram.progress(
            1 - psutil.virtual_memory().available / psutil.virtual_memory().total
        )
        cpu.progress(psutil.cpu_percent() / 100)
        time.sleep(0.5)


