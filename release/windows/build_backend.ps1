Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./build
Remove-Item -Recurse -Force -ErrorAction SilentlyContinue ./dist

# Creating the wheel
python setup.py sdist bdist_wheel
# Make sure you include the required extra packages and always use the stable or very-stable options!
pip install "dist/peptdeep-1.2.1-py3-none-any.whl[stable]"

# Creating the stand-alone pyinstaller folder
pip install pyinstaller
pyinstaller release/pyinstaller/peptdeep.spec -y


echo ls dist:
ls dist

echo ls release/windows/dist:
ls release/windows/dist

# for some reason, the installer builder expects the files here
mv dist/* release/windows/dist
mkdir release/windows/dist/peptdeep
mv release/windows/peptdeep.exe release/windows/dist/peptdeep
