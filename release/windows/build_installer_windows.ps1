# Build the installer for Windows.
# This script must be run from the root of the repository.
# Prerequisites: wheel has been build, e.g. using build_wheel.sh

Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build_pyinstaller
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist_pyinstaller


$WHL_NAME = (Get-ChildItem -Path "dist" -Filter "*.whl").Name
pip install "dist/$WHL_NAME[stable,gui-stable]"

# Creating the stand-alone pyinstaller folder
pyinstaller release/pyinstaller/peptdeep.spec  --distpath dist_pyinstaller --workpath build_pyinstaller -y
