#!/bin/bash

PACKAGE_NAME="peptdeep"
PACKAGE_VERSION="1.2.1"

ARCH=$(uname -m)
if [ "$ARCH" == "x86_64" ]; then
  ARCH="x64"
fi
echo "ARCH=${ARCH}" >> $GITHUB_ENV

KERNEL=$(uname -s | tr '[:upper:]' '[:lower:]')

BUILD_NAME="${PACKAGE_NAME}-${PACKAGE_VERSION}-${KERNEL}-${ARCH}"


# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data

# Wrapping the pyinstaller folder in a .pkg package
mkdir -p dist/peptdeep/Contents/Resources
cp release/logos/alpha_logo.icns dist/peptdeep/Contents/Resources
mv dist/peptdeep_gui dist/peptdeep/Contents/MacOS
cp release/macos/Info.plist dist/peptdeep/Contents
cp release/macos/peptdeep_terminal dist/peptdeep/Contents/MacOS
cp LICENSE.txt dist/peptdeep/Contents/Resources/LICENSE.txt
cp release/logos/alpha_logo.png dist/peptdeep/Contents/Resources/alpha_logo.png
chmod 777 release/macos/scripts/*

pkgbuild --root dist/${PACKAGE_NAME} --identifier de.mpg.biochem.${PACKAGE_NAME}.app --version 1.2.1 --install-location /Applications/${PACKAGE_NAME}.app --scripts release/macos/scripts peptdeep.pkg
productbuild --distribution release/macos/distribution.xml --resources Resources --package-path peptdeep.pkg dist/${BUILD_NAME}.pkg
