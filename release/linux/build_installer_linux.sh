#!/bin/bash
set -e -u

# Build the installer for Linux.
# This script must be run from the root of the repository.
# Prerequisites: wheel has been build, e.g. using build_wheel.sh

CPU_OR_GPU=${1:-CPU}

rm -rf dist_pyinstaller build_pyinstaller

WHL_NAME=$(cd dist && ls ./*.whl && cd ..)
pip install "dist/${WHL_NAME}[stable,gui-stable]"

if [ "${CPU_OR_GPU}" != "GPU" ]; then
  # TODO this will install the latest torch, not the version defined in stable
    pip install torch -U --extra-index-url https://download.pytorch.org/whl/cpu
fi

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/peptdeep.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
