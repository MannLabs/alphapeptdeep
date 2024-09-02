# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist

# Creating the wheel
python setup.py sdist bdist_wheel
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/peptdeep.spec -y
