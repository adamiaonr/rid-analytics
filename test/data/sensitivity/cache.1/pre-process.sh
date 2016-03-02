#!/bin/bash

usage () {
    # echo "usage: ./pre-process --scn-file <path-to-scn-file>"
    # echo "  --scn-file <path-to-scn-file>               - path to (some) .scn file"
    echo "usage: ./pre-process --scn-file <path-to-scn-file>"
    echo "  --scn-file <path-to-scn-file>               - path to .scn file"
    echo "  --help                                      - prints this menu"
}

if [[ $# -eq 0 ]]; then
    echo "ERROR: no arguments supplied"
    usage
    exit
fi

SCN_FILE=""
GRAPH_CASE="NULL"

while [[ "$1" != "" ]]; do
    
    case $1 in

        --scn-file )                    shift
                                        SCN_FILE=$1
                                        ;;
        --help )                        usage
                                        exit
                                        ;;
        * )                             usage
                                        exit
        esac
        shift

done

# get the possible cache values by first checking the tier depth on some
# .scn file (e.g. hfp.ha.scn)
tier_depth=$(cat $SCN_FILE | grep "tier_depth")
# ... and isolate it using sed to get rid of the 
# 'tier_depth=' part
tier_depth=$(echo $tier_depth | sed -e 's#.*=\(\)#\1#')

# insert a dummy line for each possible outcome (essential for the stack'd 
# graph) and <fp@>,<a@> pair

# get the possible fp probability values
fp_probs_str=$(cat $SCN_FILE | grep "fp_prob")
# ... and isolate the fp probabilities using (1) sed to get rid of the 
# 'fp_prob=' part
fp_probs_str=$(echo $fp_probs_str | sed -e 's#.*=\(\)#\1#')
# ... and the IFS (internal field separator) to separate each value into a 
# field (we save the old value of IFS in OFS, and restore it later)
OFS=$IFS
IFS=',' read -r -a fp_probs <<< "$fp_probs_str"

# do the same with alphas
alphas_str=$(cat $SCN_FILE | grep "alpha")
alphas_str=$(echo $alphas_str | sed -e 's#.*=\(\)#\1#')
IFS=',' read -r -a alphas <<< "$alphas_str"

IFS=$OFS

declare -a penalties=("feedback" "fallback")
declare -a outcomes=("correct dest." "fp detect. > orig. server" "dropped" "wrong dest. > orig. server")

for fp_prob in "${fp_probs[@]}"
do
    for alpha in "${alphas[@]}"
    do  
        for outcome in "${outcomes[@]}"
        do
            for penalty in "${penalties[@]}"
            do
                
                # add dummy values to outcomes file for each outcome, in order 
                # to make it balanced for the stack'd chart
                i=$((tier_depth))
                while [ $i -gt 0 ]
                do
                    echo "$fp_prob,$alpha,$penalty,$((i)),$outcome,0.00000000E+00,$((i - 1))" >> outcomes.csv
                    i=$((i - 1))
                done
                
            done
        done
    done
done

exit
