# copy from alphapept.gui
import os
import datetime
import yaml
import streamlit as st
from multiprocessing import Process
import psutil
import time
import pandas as pd
from typing import Callable, Union, Tuple
import pathlib

def get_posix(_path:str):
    return pathlib.PureWindowsPath(_path).as_posix()

# @st.cache
def files_in_pandas(files:list) -> pd.DataFrame:
    """Reads a folder and returns a pandas dataframe containing the files and additional information.
    Args:
        folder (str): Path to folder.

    Returns:
        pd.DataFrame: PandasDataFrame.
    """
    ctimes = [os.path.getctime(_) for _ in files]
    created = [datetime.datetime.fromtimestamp(_).strftime("%Y-%m-%d %H:%M:%S") for _ in ctimes]
    sizes = [os.path.getsize(_) / 1024 ** 2 for _ in files]
    df = pd.DataFrame(files, columns=["File Path"])
    df["Created Time"] = created
    df["File Size (Mb)"] = sizes

    return df

def file_type_selectbox( 
    ui_label:str,
    st_key:str, 
    default_type:str,
    monitor_files:list,
    choices:list,
    index=0,
)->str:
    def on_type_change():
        if len(monitor_files)>0:
            st.warning("Please clear all files before changing the file type")
            st.session_state[st_key] = default_type

    return st.selectbox(
        label=ui_label,
        options=choices, index=index,
        key=st_key,
        on_change=on_type_change
    )

def update_input_paths(file_list:list):
    _list = [
        _ for _ in file_list
        if os.path.isfile(_)
    ]
    file_list.clear()
    file_list.extend(_list)

def select_files(
    file_list:list, 
    file_exts:list, 
    ui_label="File"
):
    if isinstance(file_exts, str):
        file_exts = [file_exts.lower()]
    else:
        file_exts = [ext.lower() for ext in file_exts]
    st.write('##### ' + ui_label)
    path = st.text_input(label='Input a file or a folder path', key=ui_label+'text_input')
    path = get_posix(path)
    col1, col2, col3 = st.columns([0.5,0.5,2])
    with col1:
        add = st.button(label='Add', key=ui_label+"Add")
    with col2:
        remove = st.button(label='Remove', key=ui_label+"Remove")
    with col3:
        clear = st.button(label='Clear all files', key=ui_label+"Clear")
    if add is True:
        if os.path.isdir(path):
            for _file in os.listdir(path):
                if any(_file.lower().endswith(file_ext) for file_ext in file_exts):
                    file_list.append(os.path.join(path, _file))
        elif path not in file_list:
            file_list.append(path)
    if remove is True:
        if path in file_list:
            file_list.remove(path)
    if clear is True:
        file_list.clear()
    update_input_paths(file_list)
    st.write('##### Selected files')
    st.dataframe(files_in_pandas(file_list))

def escape_markdown(text: str) -> str:
    """Helper function to escape markdown in text.

    Args:
        text (str): Input text.

    Returns:
        str: Converted text to be used in markdown.
    """
    MD_SPECIAL_CHARS = "\`*_{}[]()#+-.!"
    for char in MD_SPECIAL_CHARS:
        text = text.replace(char, "\\" + char)
    return text


def markdown_link(description: str, link: str):
    """Creates a markdown compatible link.

    Args:
        description (str): Description.
        link (str): Target URL.
    """
    _ = f"[{description}]({link})"
    st.markdown(_, unsafe_allow_html=True)


def files_in_folder(folder: str, ending: str, sort: str = "name") -> list:
    """Reads a folder and returns all files that have this ending. Sorts the files by name or creation date.

    Args:
        folder (str): Path to folder.
        ending (str): Ending.
        sort (str, optional): How files should be sorted. Defaults to 'name'.

    Raises:
        NotImplementedError: If a sorting mode is called that is not implemented.

    Returns:
        list: List of files.
    """
    files = [_ for _ in os.listdir(folder) if _.endswith(ending)]

    if sort == "name":
        files.sort()
    elif sort == "date":
        files.sort(key=lambda x: os.path.getctime(os.path.join(folder, x)))
    else:
        raise NotImplementedError

    files = files[::-1]

    return files

def files_in_folder_pandas(folder: str, file_type:str=None) -> pd.DataFrame:
    """Reads a folder and returns a pandas dataframe containing the files and additional information.
    Args:
        folder (str): Path to folder.

    Returns:
        pd.DataFrame: PandasDataFrame.
    """
    if file_type is None:
        files = os.listdir(folder)
    else:
        file_type = file_type.lower()
        files = [
            file for file in os.listdir(folder) 
            if file.lower().endswith(f".{file_type}") or file.lower() == file_type
        ]
    created = [time.ctime(os.path.getctime(os.path.join(folder, _))) for _ in files]
    sizes = [os.path.getsize(os.path.join(folder, _)) / 1024 ** 2 for _ in files]
    df = pd.DataFrame(files, columns=["File"])
    df["Created"] = created
    df["Filesize (Mb)"] = sizes

    return df


def read_log(log_path: str):
    """Reads logfile and removes lines with __.
    Lines with __ are used to indicate progress for the AlphaPept GUI.
    Args:
        log_path (str): Path to the logile.
    """
    if os.path.isfile(log_path):
        with st.beta_expander("Run log"):
            with st.spinner("Parsing file"):
                with open(log_path, "r") as logfile:
                    lines = logfile.readlines()
                    lines = [_ for _ in lines if "__" not in _]
                    st.code("".join(lines))


def start_process(
    target: Callable,
    process_file: str,
    args: Union[list, None] = None,
    verbose: bool = True,
):
    """Function to initiate a process. It will launch the process and save the process id to a yaml file.

    Args:
        target (Callable): Target function for the process.
        process_file (str): Path to the yaml file where the process information will be stored.
        args (Union[list, None], optional): Additional arguments for the process. Defaults to None.
        verbose (bool, optional): Flag to show a stramlit message. Defaults to True.
    """
    process = {}
    now = datetime.datetime.now()
    process["created"] = now
    if args:
        p = Process(target=target, args=args)
    else:
        p = Process(target=target)
    p.start()
    process["pid"] = p.pid

    if verbose:
        st.success(f"Started process PID {p.pid} at {now}")

    with open(process_file, "w") as file:
        yaml.dump(process, file, sort_keys=False)


def check_process(
    process_path: str,
) -> Tuple[
    bool, Union[str, None], Union[str, None], Union[str, None], bool
]:
    """Function to check the status of a process.
    Reads the process file from the yaml and checks the process id.

    Args:
        process_path (str): Path to the process file.

    Returns:
        bool: Flag if process exists.
        Union ([str, None]): Process id if process exists, else None.
        Union ([str, None]): Process name if process exists, else None.
        Union ([str, None]): Process status if process exists, else None.
        bool ([type]): Flag if process was initialized.
    """
    if os.path.isfile(process_path):
        with open(process_path, "r") as process_file:
            process = yaml.load(process_file, Loader=yaml.FullLoader)
        last_pid = process["pid"]

        if "init" in process:
            p_init = process["init"]
        else:
            p_init = False

        if psutil.pid_exists(last_pid):
            p_ = psutil.Process(last_pid)
            with p_.oneshot():
                p_name = p_.name()
                status = p_.status()
            return True, last_pid, p_name, status, p_init

    return False, None, None, None, False


def init_process(process_path: str, **kwargs: dict):
    """Waits until a process file is created and then writes an init flag to the file

    Args:
        process_path (str): Path to process yaml.
    """
    while True:
        if os.path.isfile(process_path):
            with open(process_path, "r") as process_file:
                p = yaml.load(process_file, Loader=yaml.FullLoader)
            p["init"] = True
            for _ in kwargs:
                p[_] = kwargs[_]
            with open(process_path, "w") as file:
                yaml.dump(p, file, sort_keys=False)
            break
        else:
            time.sleep(1)