# navigate to the root directory

# TODO: add this to bumpversion
# OUT OF SCOPE: moving to pyproject.toml / unifying these scripts

# Creating the wheel
python setup.py sdist bdist_wheel
pip install "dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/peptdeep.spec -y
