#!/bin/bash

# Build the install package for Linux.
# This script is intended to be run from the root of the repository after running build_installer_linux.sh

PACKAGE_NAME="peptdeep"
PACKAGE_VERSION="1.2.1"
# ARCH is taken from environment variable

KERNEL=$(uname -s | tr '[:upper:]' '[:lower:]')

BUILD_NAME="${PACKAGE_NAME}-${PACKAGE_VERSION}-${KERNEL}-${ARCH}"

# If needed, include additional source such as e.g.:
# cp ../../peptdeep/data/*.fasta dist/peptdeep/data
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mv dist/peptdeep dist/peptdeep.tmp
mkdir -p dist/$BUILD_NAME/usr/local/bin
mv dist/peptdeep.tmp dist/$BUILD_NAME/usr/local/bin/peptdeep
mkdir dist/$BUILD_NAME/DEBIAN
cp release/linux/control dist/$BUILD_NAME/DEBIAN
dpkg-deb --build --root-owner-group dist/$BUILD_NAME/
