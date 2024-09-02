# Build the install package for Windows.
# This script must be run from the root of the repository after running build_installer_windows.ps1


# Apparently, ISCC uses the directory of the .iss input file as the working directory.
mv dist/* release/windows/dist
mkdir release/windows/dist/peptdeep
mv release/windows/dist/peptdeep.exe release/windows/dist/peptdeep

# Wrapping the pyinstaller folder in a .exe package
&  "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" .\release\windows\peptdeep_innoinstaller.iss
