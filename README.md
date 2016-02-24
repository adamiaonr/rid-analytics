# rid-analytics
Tool to perform simple analytical evaluations on a network using RIDs (i.e. Bloom Filters) for packet forwarding.

# Usage

Test rid-analytics quickly by running a simple example (tested on OS X El Capitan 10.11.2, 
Darwin Kernel Version 15.2.0). hopefully, the complaints from the command line will be enough to 
tell you if you need to install anything :P

## Example scenario

We will be testing the scenario depicted below

![](https://www.dropbox.com/s/v2vcngxt2t1gurc/example.png?raw=1)

Here we consider a requester R<sub>1</sub> requesting name N from the network. Content reachable via name N can be stored by producers C<sub>1,2 or 3</sub> at different points in the network.

The network is hierarchically organized in tiers. Endpoints (e.g. R<sub>1</sub> and C<sub>1,2 or 3</sub>) connect to the network via the lowest tier (t=1). The request is first looked up in routers @t1. If no positive matches happen (either false or true), the request is sent up to tier 2 and so on.

We now provide a few more characteristics for our scenario:

* 3 tiers
* 4<sup>(4-t)</sup> domains per tier
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
