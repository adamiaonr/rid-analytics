# rid-analytics
Tool to perform simple analytical evaluations on a network using RIDs (i.e. Bloom Filters) for packet forwarding.

# usage

test rid-analytics quickly by running the following sequence of commands (tested on Ubuntu 14.04.3 LTS, Linux 3.13.0-43-generic):

`$ git clone https://github.com/adamiaonr/rid-analytics.git`

`$ cd rid-analytics`

`$ make`

`$ ./rid-analytics --scn-file test/configs/sample.scn --fdist-dir test/configs/fdist --graphviz-file test/graphs/test --data-file test/data/base; dot -Tps test/graphs/test.dot -o test/graphs/test.pdf; evince test/graphs/test.pdf`

`$ cd scripts; python plot.py ../test/data/base.csv ../test/graphs/base.png`
