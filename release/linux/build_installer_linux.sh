#!/bin/bash
set -e -u

# Build the installer for Linux.
# This script must be run from the root of the repository.

CPU_OR_GPU=${1:-CPU}

rm -rf dist build *.egg-info
rm -rf dist_pyinstaller build_pyinstaller

# Creating the wheel
python -m build

# Setting up the local package
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/peptdeep-1.3.1-py3-none-any.whl[stable, gui-stable]"

if [ "${CPU_OR_GPU}" != "GPU" ]; then
    pip install torch -U --extra-index-url https://download.pytorch.org/whl/cpu
fi

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/peptdeep.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
