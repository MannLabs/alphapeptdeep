#!/bin/bash
set -e -u

# Build the wheel for all os (on a linux runner).
# This script must be run from the root of the repository.
# If you want OS-specific wheels, add the respective scripts to the OS-specific folders,
# the alphashared workflow will use those then:
# e.g. release/linux/build_wheel_linux.sh

rm -rf dist ./*.egg-info

python -m build
