# rid-analytics
Tool to perform simple analytical evaluations on a network using RIDs (i.e. Bloom Filters) for packet forwarding.

# Usage

Test rid-analytics quickly by running a simple example (tested on OS X El Capitan 10.11.2, 
Darwin Kernel Version 15.2.0). hopefully, the complaints from the command line will be enough to 
tell you if you need to install anything :P

## Example scenario

We will be testing the scenario depicted [here (slide 7)](https://www.dropbox.com/s/l466xf54fm50ax3/2016-02-14.pdf?dl=0):

* 3 tiers
* 4<sup>(4-t)</sup> domains per tier. The lowest tier, i.e. that which is closest to endpoints (content requesters and producers) has t = 1.
* We consider 4 FP rate distributions per tier:
	* H(igh) FP: All tiers 10<sup>-1</sup>
	* L(ow) FP: All tiers 10<sup>-6</sup>
	* I(ncreasing) FP (from t 0 to 3): {10<sup>-6</sup>, 10<sup>-3</sup>, 10<sup>-1</sup>}
	* D(ecreasing) FP: {10<sup>-1</sup>, 10<sup>-3</sup>, 10<sup>-6</sup>}
* 2 ALPHA distributions per tier:
	* H(igh)A: 10<sup>-1</sup>
	* L(ow)A: 10<sup>-3</sup>
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
