#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=alphadeep.pkg
if test -f "$FILE"; then
  rm alphadeep.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n alphadeepinstaller python=3.8 -y
conda activate alphadeepinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/alphadeep-0.0.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller==4.2
pyinstaller ../pyinstaller/alphadeep.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../alphadeep/data/*.fasta dist/alphadeep/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/alphadeep/Contents/Resources
cp ../logos/alpha_logo.icns dist/alphadeep/Contents/Resources
mv dist/alphadeep_gui dist/alphadeep/Contents/MacOS
cp Info.plist dist/alphadeep/Contents
cp alphadeep_terminal dist/alphadeep/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/alphadeep --identifier de.mpg.biochem.alphadeep.app --version 0.0.1 --install-location /Applications/alphadeep.app --scripts scripts alphadeep.pkg
productbuild --distribution distribution.xml --resources Resources --package-path alphadeep.pkg dist/alphadeep_gui_installer_macos.pkg
