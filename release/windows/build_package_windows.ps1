# Build the install package for Windows.
# This script must be run from the root of the repository after running build_installer_windows.ps1


# Wrapping the pyinstaller folder in a .exe package
&  "C:\Program Files (x86)\Inno Setup 6\ISCC.exe" .\release\windows\peptdeep_innoinstaller.iss
