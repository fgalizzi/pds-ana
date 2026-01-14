#!/usr/bin/env bash
# Print usage if no argument is provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 <config_file>"
    exit 1
fi

CONF_FILE=$1
echo -e "\n \n Analysing $CONF_FILE  \n \n"
eval 'root -e ".L ./VGainScans_ana.cpp+" -e "auto a = cla()" -e "VGainScans_ana(a, \"$CONF_FILE\")" -e ".q"'
