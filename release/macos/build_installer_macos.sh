#!/bin/bash
set -e -u

# Build the installer for MacOS.
# This script must be run from the root of the repository.

rm -rf dist
rm -rf build

# Creating the wheel
python setup.py sdist bdist_wheel
pip install "dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/peptdeep.spec --distpath dist_pyinstaller --workpath build_pyinstaller -y
