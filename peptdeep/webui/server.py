from concurrent.futures import process
from alphabase.yaml_utils import load_yaml

import time
import os
import shutil
import multiprocessing as mp

from peptdeep.pipeline_api import (
    generate_library, 
    transfer_learn, 
    rescore
)
from peptdeep.settings import global_settings, update_settings
from peptdeep.utils import logging

def update_global_settings(yaml_file):
    _settings = load_yaml(yaml_file)
    update_settings(global_settings, _settings)

def get_yamls(folder):
    ymls = []
    for file in os.listdir(folder):
        if file.endswith('.yaml'):
            ymls.append(os.path.join(folder, file))
    ymls.sort(key=lambda x: os.path.getmtime(x))
    return ymls

def _create_dir(dir):
    if not os.path.isdir(dir):
        os.makedirs(dir)

home_folder = os.path.expanduser(
    global_settings['PEPTDEEP_HOME']
)

queue_folder = f'{home_folder}/tasks/queue'
done_folder = f'{home_folder}/tasks/done'
failed_folder = f'{home_folder}/tasks/failed'

_create_dir(queue_folder)
_create_dir(done_folder)
_create_dir(failed_folder)

def serve():
    files = []
    echo_waiting = True
    while True:
        files = get_yamls(queue_folder)
        if len(files) > 0:
            yaml_file = files.pop(0)
            print(f"[PeptDeep] Starting a new job '{yaml_file}'...")

            try:
                update_global_settings(yaml_file)
                if global_settings['task_type'] == 'train':
                    print("[PeptDeep] Transfer learning ... ")
                    transfer_learn(global_settings)
                elif global_settings['task_type'] == 'library':
                    print("[PeptDeep] Predicting library ... ")
                    generate_library(global_settings)
                elif global_settings['task_type'] == 'rescore':
                    print("[PeptDeep] Rescoring PSMs ... ")
                    rescore(global_settings)
                else:
                    logging.warning(f"[PeptDeep] Unknown task type `{global_settings['task_type']}`, skip ... ")
                    continue
                shutil.move(
                    yaml_file, 
                    os.path.join(done_folder, os.path.basename(yaml_file))
                )
            except KeyboardInterrupt as e:
                raise e
            except Exception:
                shutil.move(
                    yaml_file, 
                    os.path.join(failed_folder, os.path.basename(yaml_file))
                )
            echo_waiting=True
        else:
            if echo_waiting:
                print("*********************************")
                print("[PeptDeep] Waiting for tasks ... ")
                print("*********************************")
                echo_waiting=False

        time.sleep(3)

class PeptDeepServer:
    def __init__(self):
        self.process:mp.Process = None

    def start(self):
        if self.process is None:
            self.process = mp.Process(target=serve)
            self.process.start()

    def terminate(self):
        if self.process is not None:
            self.process.terminate()
            self.process.kill()
            self.process = None
    
_server = PeptDeepServer()

if __name__ == '__main__':
    serve()
