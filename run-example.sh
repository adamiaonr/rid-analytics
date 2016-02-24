#!/bin/bash

usage () {
    echo "usage: bash run-example.sh --config-dir <path-to-config-dir> --data-dir <path-to-data-dir> --graph-dir <path-to-graph-dir>"
    echo "       bash run-example.sh --clean --config-dir <path-to-config-dir> --data-dir <path-to-data-dir> --graph-dir <path-to-graph-dir>"
    echo "  --config-dir <path-to-config-dir>           - path to config dir (of .scn files)"
    echo "  --data-dir <path-to-config-dir>             - path to data dir (to dump .csv and .dot files)"
    echo "  --graph-dir <path-to-graph-dir>             - path to graph dir (to save .pdf and .png files)"
    echo "  --clean                                     - just run make clean, erase all data/charts and get out"
    echo "  --help                                      - prints this menu"
}

if [[ $# -eq 0 ]]; then
    echo "ERROR: no arguments supplied"
    usage
    exit
fi

HOME_DIR=$(pwd)
CONFIG_DIR=""
DATA_DIR=""
SCRIPT_DIR="scripts"
GRAPH_DIR=""
CLEAN=0

while [[ "$1" != "" ]]; do
    
    case $1 in

        --config-dir )                  shift
                                        CONFIG_DIR=$1
                                        ;;
        --data-dir )                    shift
                                        DATA_DIR=$1
                                        ;;
        --graph-dir )                   shift
                                        GRAPH_DIR=$1
                                        ;;
        --clean )                       CLEAN=1
                                        ;;
        --help )                        usage
                                        exit
                                        ;;
        * )                             usage
                                        exit
        esac
        shift

done

declare -a fps=("hfp" "lfp" "ifp" "dfp")
declare -a alphas=("ha" "la")

# get the possible cache values by first checking the tier depth on some
# .scn file (e.g. hfp.ha.scn)
tier_depth=$(cat $CONFIG_DIR/${fps[0]}.${alphas[0]}.scn | grep "tier_depth")
# ... and isolate it using sed to get rid of the 
# 'tier_depth=' part
tier_depth=$(echo $tier_depth | sed -e 's#.*=\(\)#\1#')

# erase previous results
i=$((tier_depth))
while [ $i -gt 0 ]
do
    rm -rf $DATA_DIR/cache.$i/*.csv 
    rm -rf $DATA_DIR/cache.$i/*.dot
    rm -rf $GRAPH_DIR/*.png
    rm -rf $GRAPH_DIR/*.pdf

    i=$((i - 1))

done

if [ $CLEAN -eq 1 ]
then
    make clean
    exit
fi

for fp in "${fps[@]}"
do
    for alpha in "${alphas[@]}"
    do
        i=$((tier_depth))

        while [ $i -gt 0 ]
        do 
            ./example --scn-file $CONFIG_DIR/$fp.$alpha.scn --data-dir $DATA_DIR/cache.$i --cache-tier $((i - 1)) --test-title $fp.$alpha
            i=$((i - 1))
        done
    done
done

i=$((tier_depth))
while [ $i -gt 0 ]
do
    # pre-process the data for cache.i
    cd $DATA_DIR/cache.$i
    bash pre-process.sh --graph-dir ../../../../$GRAPH_DIR
    cd $HOME_DIR

    # run the chart scripts for it
    cd $SCRIPT_DIR
    python plot-stackd.py ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    python plot-bar.py ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    cd $HOME_DIR

    # open the .png
    open $GRAPH_DIR/stackd.cache.$i.png
    open $GRAPH_DIR/bar.cache.$i.png

    i=$((i - 1))

done

exit
