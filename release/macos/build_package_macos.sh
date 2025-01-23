#!/bin/bash
set -e -u

# Build the install package for MacOS.
# This script must be run from the root of the repository after running build_installer_macos.sh

PACKAGE_NAME=peptdeep
# BUILD_NAME is taken from environment variables, e.g. peptdeep-1.3.2-dev0-macos-darwin-arm64 or peptdeep-1.3.2-dev0-macos-darwin-x64
rm -rf ${BUILD_NAME}.pkg

# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data

# Wrapping the pyinstaller folder in a .pkg package
CONTENTS_FOLDER=dist_pyinstaller/${PACKAGE_NAME}/Contents

mkdir -p ${CONTENTS_FOLDER}/Resources
cp release/logos/alpha_logo.icns ${CONTENTS_FOLDER}/Resources
mv dist_pyinstaller/peptdeep_gui ${CONTENTS_FOLDER}/MacOS
cp release/macos/Info.plist ${CONTENTS_FOLDER}
cp release/macos/peptdeep_terminal ${CONTENTS_FOLDER}/MacOS
cp LICENSE.txt ${CONTENTS_FOLDER}/Resources
cp release/logos/alpha_logo.png ${CONTENTS_FOLDER}/Resources

# link _internal folder containing the python libraries to the Frameworks folder where they are expected
# to avoid e.g. "Failed to load Python shared library '/Applications/AlphaMap.app/Contents/Frameworks/libpython3.8.dylib'"
cd ${CONTENTS_FOLDER}
ln -s ./MacOS/_internal ./Frameworks
cd -

chmod 777 release/macos/scripts/*

pkgbuild --root dist_pyinstaller/${PACKAGE_NAME} --identifier de.mpg.biochem.${PACKAGE_NAME}.app --version 1.3.2-dev0 --install-location /Applications/${PACKAGE_NAME}.app --scripts release/macos/scripts ${PACKAGE_NAME}.pkg
productbuild --distribution release/macos/distribution.xml --resources release/macos/Resources --package-path ${PACKAGE_NAME}.pkg ${BUILD_NAME}.pkg
