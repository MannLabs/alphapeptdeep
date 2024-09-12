# Build the installer for Windows.
# This script must be run from the root of the repository.

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./*.egg-info
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build_pyinstaller
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist_pyinstaller

# Creating the wheel
python setup.py sdist bdist_wheel
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/peptdeep-1.3.0-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/peptdeep.spec  --distpath dist_pyinstaller --workpath build_pyinstaller -y
