#!/bin/bash
# This script is meant to be called by the "install" step defined in
# .travis.yml. See http://docs.travis-ci.com/ for more details.
# The behavior of the script is controlled by environment variabled defined
# in the .travis.yml in the top level folder of the project.

# License: 3-clause BSD


set -e

PYTHON_VERSION=${PYTHON_VERSION:-"2.7"}
ANUGA_PARALLEL=${ANUGA_PARALLEL:-"false"}

###########################################################
# Check if mpich2 has been installed
if [ $(dpkg-query -W -f='${Status}\n' mpich2 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="mpich2"
fi

###########################################################
# Check if openmpi has been installed
if [ $(dpkg-query -W -f='${Status}\n' openmpi-bin 2>/dev/null | grep -c "ok installed") -gt 0 ];
then
  ANUGA_PARALLEL="openmpi"
fi


echo $ANUGA_PARALLEL
