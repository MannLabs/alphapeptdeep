import time
import os
import sys
import shutil
import multiprocessing as mp

from peptdeep.pipeline_api import (
    generate_library, 
    transfer_learn, 
    rescore
)
from peptdeep.settings import global_settings, load_global_settings
from peptdeep.utils import logging

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
            
            running_txt = f'{home_folder}/tasks/running.txt'
            with open(running_txt,'w') as f:
                f.write(yaml_file)

            try:
                load_global_settings(yaml_file)
                if global_settings['task_type'] == 'train':
                    print("[PeptDeep] Transfer learning ... ")
                    transfer_learn()
                elif global_settings['task_type'] == 'library':
                    print("[PeptDeep] Predicting library ... ")
                    generate_library()
                elif global_settings['task_type'] == 'rescore':
                    print("[PeptDeep] Rescoring PSMs ... ")
                    rescore()
                else:
                    logging.warning(f"[PeptDeep] Unknown task type `{global_settings['task_type']}`, skip ... ")
                    continue
                if os.path.isfile(yaml_file):
                    shutil.move(
                        yaml_file, 
                        os.path.join(done_folder, os.path.basename(yaml_file))
                    )
            except KeyboardInterrupt as e:
                with open(running_txt,'w') as f:
                    f.write("")
                raise e
            except Exception as e:
                if os.path.isfile(yaml_file):
                    shutil.move(
                        yaml_file, 
                        os.path.join(failed_folder, os.path.basename(yaml_file)),
                    )
                print(e)
            with open(running_txt,'w') as f:
                f.write("")
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
        self._process_file = os.path.join(
            global_settings['PEPTDEEP_HOME'], 
            'tasks/serve_pid.txt', 
        )

    def start(self):
        if self.process is None:
            self.process = mp.Process(target=serve)
            self.process.start()

            with open(self._process_file, 'w') as f:
                f.write(str(self.process.pid))

    def terminate(self):
        if self.process is not None:
            self.process.terminate()
            self.process.kill()
            self.process = None

        os.replace(self._process_file, self._process_file[:-3]+'prev.txt')

    def __del__(self):
        self.terminate()
    
_server = PeptDeepServer()
