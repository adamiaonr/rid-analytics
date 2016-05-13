#!/bin/bash

usage () {
    echo "usage: bash run-sanity.sh --test-case <case-name> --bf-size <size-in-bits> --fwd-table-size <nr-of-entreis> --scn-dir <path-to-scn-dir> --data-dir <path-to-data-dir> --graph-dir <path-to-graph-dir>"
    echo "       bash run-sanity.sh --clean --scn-dir <path-to-scn-dir> --data-dir <path-to-data-dir> --graph-dir <path-to-graph-dir>"
    echo "  --test-case                                 - the test case to run. e.g.: 'no-caching', 'prefix-management'"
    echo "  --bf-size <size-in-bits>                    - size of RID (Bloom filter) to use in the evaluations"
    echo "  --fwd-table-size <nr-of-entreis>            - nr. of RID fowarding entries to use in base forwarding table"
    echo "  --scn-dir <path-to-scn-dir>                 - path to .scn file dir. this script will read all .scn files and run evaluations for them."
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

# save home dir (i.e. where run-analysis.sh is being called)
HOME_DIR=$(pwd)
CACHING="caching"
PERVASIVE_CACHING="pervasive"
EDGE_CACHING="edge"

# input parameter placeholders
SCN_DIR=""
TEST_CASE=""
BF_SIZE=160
FWD_TABLE_SIZE=100000000
SCN_DIR=""
SCN_DIR_ALT="test/data/sanity-check/default"
DATA_DIR=""
DATA_SUB_DIR=(5 10)
EXPAND_FACTORS=(1 2 5 10)
GRAPH_DIR=""
SCRIPT_DIR="scripts"

# we don't clean by default (obviously)
CLEAN=0

while [[ "$1" != "" ]]; do
    
    case $1 in

        --test-case )                   shift
                                        TEST_CASE=$1
                                        ;;
        --bf-size )                     shift
                                        BF_SIZE=$1
                                        ;;
        --fwd-table-size )              shift
                                        FWD_TABLE_SIZE=$1
                                        ;;
        --scn-dir )                     shift
                                        SCN_DIR=$1
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

if [ $CLEAN -eq 1 ]
then
    make clean
    exit
fi

if [ "${TEST_CASE}" = 'no-caching' ]; then

    for sub_dir in ${DATA_SUB_DIR[@]} 
    do
        for file in $SCN_DIR/r$sub_dir/*.scn
        do
            echo "processing $file"

            # strip the '.scn' extension from the file name 
            file_name=$(basename "$file")
            file_name="${file_name%.*}"-$BF_SIZE

            # run the sanity check for each test case
            ./sanity-check --bf-size $BF_SIZE --request-size $sub_dir --fwd-table-size $FWD_TABLE_SIZE --scn-file $file --data-dir $DATA_DIR/r$sub_dir --output-file $file_name -v
        done
    done

    cd $SCRIPT_DIR
    python plot-outcomes.py $BF_SIZE $DATA_DIR $GRAPH_DIR no-caching-$BF_SIZE

    cd $HOME_DIR
    open $GRAPH_DIR/no-caching-$BF_SIZE.pdf
fi

i=1
if [ "${TEST_CASE}" = 'prefix-management' ]; then

    while [ $i -le 10 ]
    do
        file=$SCN_DIR/r10/pm-5-10-th$i.scn
        echo "processing $file"

        # strip the '.scn' extension from the file name 
        file_name=$(basename "$file")

        for exp in "${EXPAND_FACTORS[@]}"
        do
            output_file_name="${file_name%.*}"-$BF_SIZE-$i-$exp
            # run the sanity check for each test case
            ./sanity-check --f-min $i --expand-factor $exp --bf-size $BF_SIZE --request-size 10 --fwd-table-size $FWD_TABLE_SIZE --scn-file $file --data-dir $DATA_DIR/r10 --output-file $output_file_name -v
        done

        i=$((i + 1))

    done

    cd $SCRIPT_DIR
    python plot-prefix-management.py $BF_SIZE $DATA_DIR $GRAPH_DIR prefix-management-5-10-$BF_SIZE

    cd $HOME_DIR
    open $GRAPH_DIR/prefix-management-5-10-$BF_SIZE.pdf
fi

if [ "${TEST_CASE}" = 'on-path' ]; then

    for sub_dir in ${DATA_SUB_DIR[@]} 
    do

        for scn_file in $SCN_DIR/$CACHING/$PERVASIVE_CACHING/on-path/r$sub_dir/*.scn
        do
            echo "processing $scn_file"

            # strip the '.scn' extension from the file name 
            output_file_name=$(basename "$scn_file")
            output_file_name="${output_file_name%.*}"-$BF_SIZE

            # output data dir
            output_dir=$DATA_DIR/$CACHING/$PERVASIVE_CACHING/on-path/r$sub_dir

            # run the sanity check for each test case
            ./sanity-check --bf-size $BF_SIZE --request-size $sub_dir --fwd-table-size $FWD_TABLE_SIZE --scn-file $scn_file --data-dir $output_dir --output-file $output_file_name -v
        done

    done
fi

if [ "${TEST_CASE}" = 'nearest-replica' ]; then

    for sub_dir in ${DATA_SUB_DIR[@]} 
    do

        for scn_file in $SCN_DIR/$CACHING/$PERVASIVE_CACHING/nearest-replica/r$sub_dir/*.scn
        do
            echo "processing $scn_file"

            # strip the '.scn' extension from the file name 
            output_file_name=$(basename "$scn_file")
            output_file_name="${output_file_name%.*}"-$BF_SIZE

            # output data dir
            output_dir=$DATA_DIR/$CACHING/$PERVASIVE_CACHING/nearest-replica/r$sub_dir

            # run the sanity check for each test case
            ./sanity-check --bf-size $BF_SIZE --request-size $sub_dir --fwd-table-size $FWD_TABLE_SIZE --scn-file $scn_file --data-dir $output_dir --output-file $output_file_name -v
        done

    done
fi

if [ "${TEST_CASE}" = 'coop' ]; then

    for sub_dir in ${DATA_SUB_DIR[@]} 
    do

        for scn_file in $SCN_DIR/$CACHING/$EDGE_CACHING/coop/r$sub_dir/*.scn
        do
            echo "processing $scn_file"

            # strip the '.scn' extension from the file name 
            output_file_name=$(basename "$scn_file")
            output_file_name="${output_file_name%.*}"-$BF_SIZE

            # output data dir
            output_dir=$DATA_DIR/$CACHING/$EDGE_CACHING/coop/r$sub_dir

            # run the sanity check for each test case
            ./sanity-check --bf-size $BF_SIZE --request-size $sub_dir --fwd-table-size $FWD_TABLE_SIZE --scn-file $scn_file --data-dir $output_dir --output-file $output_file_name -v
        done

    done
fi

if [ "${TEST_CASE}" = 'offpath-coop' ]; then

    for sub_dir in ${DATA_SUB_DIR[@]} 
    do

        for scn_file in $SCN_DIR/$CACHING/$EDGE_CACHING/offpath-coop/r$sub_dir/*.scn
        do
            echo "processing $scn_file"

            # strip the '.scn' extension from the file name 
            output_file_name=$(basename "$scn_file")
            output_file_name="${output_file_name%.*}"-$BF_SIZE

            # output data dir
            output_dir=$DATA_DIR/$CACHING/$EDGE_CACHING/offpath-coop/r$sub_dir

            # run the sanity check for each test case
            ./sanity-check --bf-size $BF_SIZE --request-size $sub_dir --fwd-table-size $FWD_TABLE_SIZE --scn-file $scn_file --data-dir $output_dir --output-file $output_file_name -v
        done

    done
fi

exit