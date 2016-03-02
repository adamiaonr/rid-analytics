#!/bin/bash

parameter_change () {

    key=$1
    value=$2
    filename=$3

    echo $key
    echo $value
    echo $filename

    # let's make this complex and use sed to change the value of key 'in-place'
    sed -i '' "s/^$1.*$/$1$2/" "$3"
}

usage () {
    echo "usage: bash run-analysis.sh --scn-file <path-to-scn-file> --data-dir <path-to-data-dir> --graph-dir <path-to-graph-dir>"
    echo "       bash run-analysis.sh --clean --scn-file <path-to-scn-file> --data-dir <path-to-data-dir> --graph-dir <path-to-graph-dir>"
    echo "  --scn-file <path-to-scn-file>               - path to .scn file"
    echo "  --data-dir <path-to-config-dir>             - path to data dir (to dump .csv and .dot files)"
    echo "  --graph-dir <path-to-graph-dir>             - path to graph dir (to save .pdf and .png files)"
    echo "  --clean                                     - just run make clean, erase all data/charts and get out"
    echo "  --quick-change <parameter> <value>          - quickly change one of the test parameters. e.g. '--quick-change ha 0.5,0.5,0.5' changes the value of every 'alpha=' parameter on xfp.ha.scn files. optional and may be used multiple times."
    echo "  --help                                      - prints this menu"
}

if [[ $# -eq 0 ]]; then
    echo "ERROR: no arguments supplied"
    usage
    exit
fi

# save home dir (i.e. where run-analysis.sh is being called)
HOME_DIR=$(pwd)

# input parameter placeholders
SCN_FILE=""
SCN_FILE_ALT="test/configs/sensitivity/sensitivity.dist.scn"
DATA_DIR=""
GRAPH_DIR=""
SCRIPT_DIR="scripts"

# we don't clean by default (obviously)
CLEAN=0

# placeholders for quick parameter changes
PARAMETER_CHANGE=0
PARAMETER_KEYS=()
PARAMETER_VALUES=()

while [[ "$1" != "" ]]; do
    
    case $1 in

        --scn-file )                    shift
                                        SCN_FILE=$1
                                        ;;
        --data-dir )                    shift
                                        DATA_DIR=$1
                                        ;;
        --graph-dir )                   shift
                                        GRAPH_DIR=$1
                                        ;;
        --quick-change )                PARAMETER_CHANGE=1
                                        shift
                                        PARAMETER_KEYS+=("$1")
                                        shift
                                        PARAMETER_VALUES+=("$1")
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

i=0
if [ $PARAMETER_CHANGE -eq 1 ]
then
    for key in "${PARAMETER_KEYS[@]}"
    do
        for alpha in "${alphas[@]}"
        do
            if [ $alpha = $key ]
            then
                for fp in "${fps[@]}"
                do
                    # call parameter change of <fp>.<alpha>.scn
                    parameter_change "alpha=" ${PARAMETER_VALUES[$i]} $SCN_FILE

                done
            fi
        done

        i=$((i + 1))
    done
fi

# get the possible cache values by first checking the tier depth on some
# .scn file (e.g. hfp.ha.scn)
tier_depth=$(cat $SCN_FILE | grep "tier_depth")
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

i=$((tier_depth))

while [ $i -gt 0 ]
do 
    ./sensitivity-analysis --scn-file $SCN_FILE --data-dir $DATA_DIR/cache.$i --cache-tier $((i - 1)) --test-title cache.$i
    i=$((i - 1))
done

i=$((tier_depth))
while [ $i -gt 0 ]
do
    # pre-process the data for cache.i
    cd $DATA_DIR/cache.$i
    bash pre-process.sh --scn-file ../../../../$SCN_FILE
    cd $HOME_DIR

    # run the chart scripts for it
    cd $SCRIPT_DIR

    python plot-outcome-vs-alpha.py ../$SCN_FILE ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    python plot-outcome-vs-tier.py ../$SCN_FILE_ALT ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    python plot-outcome-vs-hop.py ../$SCN_FILE_ALT ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    # python plot-stackd.py fprob $tier_depth ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    # python plot-stackd.py alpha $tier_depth ../$DATA_DIR/cache.$i ../$GRAPH_DIR

    # python plot-box.py fprob $tier_depth ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    # python plot-box.py alpha $tier_depth ../$DATA_DIR/cache.$i ../$GRAPH_DIR
    cd $HOME_DIR

    # open the .png
    open $GRAPH_DIR/outcome-vs-alpha.fallback.cache.$i.png
    open $GRAPH_DIR/outcome-vs-tier.fallback.cache.$i.png
    open $GRAPH_DIR/outcome-vs-hop.fallback.cache.$i.png
    # open $GRAPH_DIR/stackd.fprob.fallback.cache.$i.png
    # open $GRAPH_DIR/stackd.alpha.fallback.cache.$i.png
    # open $GRAPH_DIR/box.fprob.cache.$i.png
    # open $GRAPH_DIR/box.fprob.cache.$i.png

    i=$((i - 1))

done

cd $SCRIPT_DIR
python plot-outcome-vs-dist.py ../$SCN_FILE_ALT ../$DATA_DIR ../$GRAPH_DIR
cd $HOME_DIR

open $GRAPH_DIR/outcome-vs-dist.fallback.png

exit
