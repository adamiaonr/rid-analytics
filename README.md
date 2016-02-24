# rid-analytics
Tool to perform simple analytical evaluations on a network using RIDs (i.e. Bloom Filters) for packet forwarding.

# Usage

Test rid-analytics quickly by running a simple example (tested on OS X El Capitan 10.11.2, 
Darwin Kernel Version 15.2.0). hopefully, the complaints from the command line will be enough to 
tell you if you need to install anything :P

## Example scenario

We will be testing the scenario depicted [here](https://www.dropbox.com/s/v2vcngxt2t1gurc/example.png?dl=0):

* 3 tiers
* 4^(4-t) domains per tier. The lowest tier, i.e. that which is closest to endpoints (content requesters and producers) has t = 1.
* We consider 4 FP rate distributions per tier:
	* H(igh) FP: All tiers 10^-1
	* L(ow) FP: All tiers 10^-6
	* I(ncreasing) FP (from t 0 to 3): {10^-6, 10^-3, 10^-1}
	* D(ecreasing) FP: {10^-1, 10^-3, 10^-6}
* 2 ALPHA distributions per tier:
	* H(igh)A: 10^-1
	* L(ow)A: 10^-3
* 3 cases for cache locations {# of content sources @t1,@t2,@t3}: 
	* {1,1,1}
	* {0,1,1} 
	* {0,0,1}	

## Running the example

Open Terminal and run the following list of commands:

`$ git clone https://github.com/adamiaonr/rid-analytics.git`

`$ cd rid-analytics`

`$ git checkout present`

`$ make`

`$ bash run-example.sh --clean --config-dir test/configs/example --data-dir test/data/example --graph-dir test/graphs/example`

## Output

(...) 
