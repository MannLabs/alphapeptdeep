# Build the installer for Windows.
# This script is intended to be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist

# temporary fix for the issue with windows installer getting stuck
pip install numpy==2.0.1 altair==5.4.0 packaging==24.1 rich==13.7.1 idna==3.7 importlib_metadata==8.4.0 zipp==3.20.0

# Creating the wheel
python setup.py sdist bdist_wheel
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/peptdeep.spec -y

# for some reason, the installer builder expects the files here
mv dist/* release/windows/dist
mkdir release/windows/dist/peptdeep
mv release/windows/dist/peptdeep.exe release/windows/dist/peptdeep
