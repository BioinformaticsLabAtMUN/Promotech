#!/bin/bash 
set -e 

export bTSSfinder_Data=$4
# echo "PATH: " $PWD
$1 -i $2 -o $3 -h 2 -c 1 -t $5