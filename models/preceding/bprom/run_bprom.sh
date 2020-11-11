#!/bin/bash 

set -e 
# SETTING DEPENDENCIES PATH
# export TSS_DATA="$(pwd)/bprom_data"
export TSS_DATA="$5"


# RUNNING BPROM

# echo "RUNNING BPROM IN FILE: " $1 
# echo "SAVING IN " $2 
# echo "TSS_DATA PATH: " $TSS_DATA 
# echo "CMD: " $4 $1 $2
# ./bprom $1 $2
$4 $1 $2

# PARSING OUPUT FILE
# grep -E '>|LDF|-10|-35' $2 | awk -F "\n" '!/^>/ {printf "\t%s", $0; n= "\n"} /^>/ {printf "\n%s", $0; n = ""} END {printf "\t%s", n}' | grep Promoter | perl -p -e 's/ +/\t/g' | perl -p -e 's/\t+/\t/g' | cut -f1,4,6,11,12,14,19,20,22 
grep -E '>|LDF|-10|-35' $2 | awk -F "\n" '!/^>/ {printf "\t%s", $0; n= "\n"} /^>/ {printf "\n%s", $0; n = ""} END {printf "\t%s", n}' | grep Promoter | perl -p -e 's/ +/\t/g' | perl -p -e 's/\t+/\t/g' | cut -f1,4,6,11,12,14,19,20,22 > $3
