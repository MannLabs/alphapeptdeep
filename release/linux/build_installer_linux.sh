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
# WARNING: this probably does not work!!!!

# Wrapping the pyinstaller folder in a .deb package
mv dist/peptdeep mv dist/peptdeep.tmp
mkdir -p dist/$BUILD_NAME/usr/local/bin
mv dist/peptdeep.tmp dist/$BUILD_NAME/usr/local/bin/peptdeep
mkdir dist/$BUILD_NAME/DEBIAN
cp release/linux/control dist/$BUILD_NAME/DEBIAN
dpkg-deb --build --root-owner-group dist/$BUILD_NAME/
