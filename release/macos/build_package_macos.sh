#!/bin/bash
set -u -e

# Build the install package for MacOS.
# This script is intended to be run from the root of the repository after running build_installer_macos.sh

PACKAGE_NAME=peptdeep
# BUILD_NAME is taken from environment variables, e.g. peptdeep-1.2.1-macos-darwin-arm64.pkg or peptdeep-1.2.1-macos-darwin-x64.pkg
rm -rf ${BUILD_NAME}.pkg

# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/${PACKAGE_NAME}/Contents/Resources
cp release/logos/alpha_logo.icns dist/${PACKAGE_NAME}/Contents/Resources
mv dist/peptdeep_gui dist/${PACKAGE_NAME}/Contents/MacOS
cp release/macos/Info.plist dist/${PACKAGE_NAME}/Contents
cp release/macos/peptdeep_terminal dist/${PACKAGE_NAME}/Contents/MacOS
cp LICENSE.txt dist/${PACKAGE_NAME}/Contents/Resources/LICENSE.txt
cp release/logos/alpha_logo.png dist/${PACKAGE_NAME}/Contents/Resources/alpha_logo.png
chmod 777 release/macos/scripts/*

pkgbuild --root dist/${PACKAGE_NAME} --identifier de.mpg.biochem.${PACKAGE_NAME}.app --version 1.2.1 --install-location /Applications/${PACKAGE_NAME}.app --scripts release/macos/scripts ${PACKAGE_NAME}.pkg
productbuild --distribution release/macos/distribution.xml --resources Resources --package-path ${PACKAGE_NAME}.pkg dist/${BUILD_NAME}.pkg
