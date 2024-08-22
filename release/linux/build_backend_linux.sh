#!/bin/bash

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

if [ "$1" == "CPU" ]; then
    pip install torch -U --extra-index-url https://download.pytorch.org/whl/cpu
fi

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/peptdeep.spec -y
