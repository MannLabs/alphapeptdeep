#!bash

# Initial cleanup
rm -rf dist
rm -rf build
FILE=peptdeep.pkg
if test -f "$FILE"; then
  rm peptdeep.pkg
fi
cd ../..
rm -rf dist
rm -rf build

# Creating a conda environment
conda create -n peptdeepinstaller python=3.8 -y
conda activate peptdeepinstaller

# Creating the wheel
python setup.py sdist bdist_wheel

# Setting up the local package
cd release/one_click_macos_gui
pip install "../../dist/peptdeep-0.4.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller ../pyinstaller/peptdeep.spec -y
conda deactivate

# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/peptdeep/Contents/Resources
cp ../logos/alpha_logo.icns dist/peptdeep/Contents/Resources
mv dist/peptdeep_gui dist/peptdeep/Contents/MacOS
cp Info.plist dist/peptdeep/Contents
cp peptdeep_terminal dist/peptdeep/Contents/MacOS
cp ../../LICENSE.txt Resources/LICENSE.txt
cp ../logos/alpha_logo.png Resources/alpha_logo.png
chmod 777 scripts/*

pkgbuild --root dist/peptdeep --identifier de.mpg.biochem.peptdeep.app --version 0.4.1 --install-location /Applications/peptdeep.app --scripts scripts peptdeep.pkg
productbuild --distribution distribution.xml --resources Resources --package-path peptdeep.pkg dist/peptdeep_gui_installer_macos.pkg
