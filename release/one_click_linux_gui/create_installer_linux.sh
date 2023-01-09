#!bash

# Initial cleanup
rm -rf dist
rm -rf build
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n peptdeep_installer python=3.10 -y
conda activate peptdeep_installer

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_linux_gui
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "../../dist/peptdeep-1.0.0-py3-none-any.whl[stable]"

if [ "$1" == "CPU" ]; then
    pip install torch -U --extra-index-url https://download.pytorch.org/whl/cpu
fi

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller ../pyinstaller/peptdeep.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mkdir -p dist/peptdeep_gui_installer_linux/usr/local/bin
mv dist/peptdeep dist/peptdeep_gui_installer_linux/usr/local/bin/peptdeep
mkdir dist/peptdeep_gui_installer_linux/DEBIAN
cp control dist/peptdeep_gui_installer_linux/DEBIAN
dpkg-deb --build --root-owner-group dist/peptdeep_gui_installer_linux/
