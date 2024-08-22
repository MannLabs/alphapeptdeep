#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n peptdeep_installer python=3.9 -y
conda activate peptdeep_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/windows
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller ../pyinstaller/peptdeep.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data

# Wrapping the pyinstaller folder in a .exe package
"C:\Program Files (x86)\Inno Setup 6\ISCC.exe" peptdeep_innoinstaller.iss
# WARNING: this assumes a static location for innosetup
