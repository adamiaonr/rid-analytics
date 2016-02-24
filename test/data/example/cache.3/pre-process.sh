#!/bin/bash

usage () {
    # echo "usage: ./pre-process --scn-file <path-to-scn-file>"
    # echo "  --scn-file <path-to-scn-file>               - path to (some) .scn file"
    echo "usage: ./pre-process --graph-dir <path-to-graph-dir> --show-graph <fp.alpha> (optional arg.)"
    echo "  --graph-dir <path-to-graph-dir>             - path to graph dir (to save .pdf and .png files)"
    echo "  --show-graph <fp.alpha>                     - show probability tree for <fp.alpha> case (optional)"
    echo "  --help                                      - prints this menu"
}

# if [[ $# -eq 0 ]]; then
#     echo "ERROR: no arguments supplied"
#     usage
#     exit
# fi

SCN_FILE=""
GRAPH_DIR=""
GRAPH_CASE="NULL"

while [[ "$1" != "" ]]; do
    
    case $1 in

        --graph-dir )                   shift
                                        GRAPH_DIR=$1
                                        ;;
        --show-graph )                  shift
                                        GRAPH_CASE=$1
                                        ;;
        --help )                        usage
                                        exit
                                        ;;
        * )                             usage
                                        exit
        esac
        shift

done

# # get the possible fp probability values
# fp_prob_str=$(cat $SCN_FILE | grep "fp_prob")
# # ... and isolate the fp probabilities using (1) sed to get rid of the 
# # 'fp_prob=' part
# fp_prob_str=$(echo $fp_prob_str | sed -e 's#.*=\(\)#\1#')
# # ... and the IFS (internal field separator) to separate each value into a 
# # field (we save the old value of IFS in OFS, and restore it later)
# OFS=$IFS
# IFS=','
# fp_prob=$fp_prob_str
# IFS=$OFS

# insert a dummy line for each possible outcome (essential for the stack'd 
# graph)
declare -a fps=("hfp" "lfp" "ifp" "dfp")
declare -a alphas=("ha" "la")
declare -a penalties=("feedback" "fallback")

declare -a outcome_str=("correct dest." "fp detect. > orig. server" "dropped" "wrong dest. > orig. server")
declare -a outcome_short_str=("cc" "rlyd" "drpd" "ic")

for fp in "${fps[@]}"
do
    for alpha in "${alphas[@]}"
    do
        index=0

        dot -Tpdf $fp.$alpha.dot -o $GRAPH_DIR/$fp.$alpha.pdf

        for outcome in "${outcome_str[@]}"
        do
            for penalty in "${penalties[@]}"
            do
                # separate latencies into files according to penalty type
                cat latencies.csv | grep "$penalty" > latencies.$penalty.csv

                # add dummy values to outcomes file for each outcome, in order 
                # to make it balanced for the stack'd chart
                echo "$fp.$alpha,$penalty,$outcome,0.00000000E+00" >> outcomes.csv

                # separate outcomes according to outcome and penalty type
                cat outcomes.csv | grep "$penalty" | grep "$outcome" > outcomes.$penalty.${outcome_short_str[$index]}.csv
            done

            index=$((index + 1))
        done
    done
done

if [ $GRAPH_CASE != "NULL" ]
then
    open $GRAPH_DIR/$GRAPH_CASE.pdf
fi

exit
