#!/bin/bash

# License: 3-clause BSD


set -e

echo "#==========================="
echo "# Determine Ubuntu version"
echo "#==========================="

VERSION_ID=$(grep -oP 'VERSION_ID="\K[\d.]+' /etc/os-release)

echo "Version $VERSION_ID"
echo " "

if [[ "$VERSION_ID" == "20.04" ]] 
then 
    cd "$(dirname "${BASH_SOURCE[0]}")";
    bash install_ubuntu_20_04.sh
fi

if [[ "$VERSION_ID" == "22.04" ]] 
then 
    cd "$(dirname "${BASH_SOURCE[0]}")";
    bash install_ubuntu_22_04.sh
fi

if [[ "$VERSION_ID" == "24.04" ]] 
then 
    cd "$(dirname "${BASH_SOURCE[0]}")";
    bash install_ubuntu_24_04.sh
fi