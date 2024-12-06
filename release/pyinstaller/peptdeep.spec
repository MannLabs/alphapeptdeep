# -*- mode: python ; coding: utf-8 -*-

import pkgutil
import os
import sys
from PyInstaller.building.build_main import Analysis, PYZ, EXE, COLLECT, BUNDLE, TOC
import PyInstaller.utils.hooks
from PyInstaller.utils.hooks import copy_metadata
import pkg_resources
import importlib.metadata
from transformers.dependency_versions_check import pkgs_to_check_at_runtime
import peptdeep
from peptdeep.utils._pyinstaller_hooks import get_peptdeep_datas


##################### User definitions
gui_name = 'peptdeep_gui'
gui_script = 'peptdeep_pyinstaller.py'
cli_name = 'peptdeep'
cli_script = 'peptdeep_cli_pyinstaller.py'
if sys.platform[:6] == "darwin":
	icon = '../logos/alpha_logo.icns'
else:
	icon = '../../peptdeep/webui/logos/peptdeep.ico'
block_cipher = None
location = os.getcwd()
project = "peptdeep"
remove_tests = True
bundle_name = "peptdeep"
#####################


datas, binaries, hidden_imports = PyInstaller.utils.hooks.collect_all(
	project,
	include_py_files=True
)

alpha_x = ['alphabase', 'alpharaw']
for alpha_package in alpha_x:
	_datas, _binaries, _hidden_imports = PyInstaller.utils.hooks.collect_all(
		alpha_package,
		include_py_files=True
	)
	datas+=_datas
	binaries+=_binaries
	hidden_imports+=_hidden_imports


hidden_imports = [h for h in hidden_imports if "__pycache__" not in h]
datas = [d for d in datas if ("__pycache__" not in d[0]) and (d[1] not in [".", "Resources", "scripts"])]

#if sys.platform[:5] == "win32":
#	base_path = os.path.dirname(sys.executable)
#	library_path = os.path.join(base_path, "Library", "bin")
#	dll_path = os.path.join(base_path, "DLLs")
#	libcrypto_dll_path = os.path.join(dll_path, "libcrypto-1_1-x64.dll")
#	libssl_dll_path = os.path.join(dll_path, "libssl-1_1-x64.dll")
#	libcrypto_lib_path = os.path.join(library_path, "libcrypto-1_1-x64.dll")
#	libssl_lib_path = os.path.join(library_path, "libssl-1_1-x64.dll")
#	if not os.path.exists(libcrypto_dll_path):
#		datas.append((libcrypto_lib_path, "."))
#	if not os.path.exists(libssl_dll_path):
#		datas.append((libssl_lib_path, "."))

datas += get_peptdeep_datas()

gui_a = Analysis(
	[gui_script],
	pathex=[location],
	binaries=binaries,
	datas=datas,
	hiddenimports=hidden_imports,
	hookspath=[],
	runtime_hooks=[],
	excludes=[h for h in hidden_imports if "datashader" in h],
	win_no_prefer_redirects=False,
	win_private_assemblies=False,
	cipher=block_cipher,
	noarchive=False
)

gui_pyz = PYZ(
	gui_a.pure,
	gui_a.zipped_data,
	cipher=block_cipher
)

if sys.platform.startswith("linux"):
	gui_exe = EXE(
		gui_pyz,
		gui_a.scripts,
		gui_a.binaries,
		gui_a.zipfiles,
		gui_a.datas,
		name=bundle_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		upx_exclude=[],
		icon=icon
	)
elif sys.platform.startswith('win32'):
	cli_a = Analysis(
		[cli_script],
		pathex=[location],
		binaries=binaries,
		datas=datas,
		hiddenimports=hidden_imports,
		hookspath=[],
		runtime_hooks=[],
		excludes=[h for h in hidden_imports if "datashader" in h],
		win_no_prefer_redirects=False,
		win_private_assemblies=False,
		cipher=block_cipher,
		noarchive=False
	)
	cli_pyz = PYZ(
		cli_a.pure,
		cli_a.zipped_data,
		cipher=block_cipher
	)
	MERGE( (gui_a, 'gui', 'gui'), (cli_a, 'cli', 'cli') )

	gui_exe = EXE(
		gui_pyz,
		gui_a.scripts,
		# a.binaries,
		gui_a.zipfiles,
		# a.datas,
		exclude_binaries=True,
		name=gui_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		icon=icon
	)
	gui_coll = COLLECT(
		gui_exe,
		gui_a.binaries,
		# a.zipfiles,
		gui_a.datas,
		strip=False,
		upx=True,
		upx_exclude=[],
		name=gui_name
	)

	cli_exe = EXE(
		cli_pyz,
		cli_a.scripts,
		# a.binaries,
		cli_a.zipfiles,
		# a.datas,
		exclude_binaries=True,
		name=cli_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		icon=icon
	)
	cli_coll = COLLECT(
		cli_exe,
		cli_a.binaries,
		# a.zipfiles,
		cli_a.datas,
		strip=False,
		upx=True,
		upx_exclude=[],
		name=cli_name
	)
else:  # macOS
	gui_exe = EXE(
		gui_pyz,
		gui_a.scripts,
		# a.binaries,
		gui_a.zipfiles,
		# a.datas,
		exclude_binaries=True,
		name=gui_name,
		debug=False,
		bootloader_ignore_signals=False,
		strip=False,
		upx=True,
		console=True,
		icon=icon
	)
	gui_coll = COLLECT(
		gui_exe,
		gui_a.binaries,
		# a.zipfiles,
		gui_a.datas,
		strip=False,
		upx=True,
		upx_exclude=[],
		name=gui_name
	)
