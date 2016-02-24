# rid-analytics
<a name="sec:intro"></a>

Tool to perform simple analytical evaluations on a network using RIDs (i.e. Bloom Filters) for packet forwarding.

# Usage
<a name="sec:usage"></a>

Test rid-analytics quickly by running a simple example (tested on OS X El Capitan 10.11.2, 
Darwin Kernel Version 15.2.0). hopefully, the complaints from the command line will be enough to 
tell you if you need to install anything :P

## Example scenario
<a name="subsec:scn"></a>

We will be testing the scenario depicted below

![](https://www.dropbox.com/s/v2vcngxt2t1gurc/example.png?raw=1)

Here we consider a requester R<sub>1</sub> requesting name N from the network. Content reachable via name N can be stored by producers C<sub>1,2 or 3</sub> at different points in the network.

The network is hierarchically organized in tiers. Endpoints (e.g. R<sub>1</sub> and C<sub>1,2 or 3</sub>) connect to the network via the lowest tier (t1 or 'tier 1'). The request is first looked up in a router @t1. If no positive matches happen (either false or true), the request is sent up to tier 2 and so on.

A few more characteristics for our scenario:

* 3 tiers
* 4<sup>(4-t)</sup> domains per tier
* We consider 4 FP rate distributions per tier:
	* H(igh) FP: All tiers 10<sup>-1</sup>
	* L(ow) FP: All tiers 10<sup>-6</sup>
	* I(ncreasing) FP (from t1 to t3): {10<sup>-6</sup>, 10<sup>-3</sup>, 10<sup>-1</sup>}
	* D(ecreasing) FP: {10<sup>-1</sup>, 10<sup>-3</sup>, 10<sup>-6</sup>}
* 2 ALPHA distributions per tier:
	* H(igh)A: 10<sup>-1</sup>
	* L(ow)A: 10<sup>-3</sup>
* 3 cases for cache locations {# of content sources @t1,@t2,@t3}: 
	* {1,1,1}
	* {0,1,1} 
	* {0,0,1}	

## Running the example
<a name="subsec:run"></a>

Open Terminal and run the following list of commands:

`$ git clone https://github.com/adamiaonr/rid-analytics.git`

`$ cd rid-analytics`

`$ git checkout present`

`$ make`

`$ bash run-example.sh --config-dir test/configs/example --data-dir test/data/example --graph-dir test/graphs/example`

## Output
<a name="subsec:output"></a>

After running `run-example` you should have some .png and pdf files on `test/graphs/example/`. Look down for a description.

### Outcomes
<a name="subsubsec:outcomes"></a>

![](https://www.dropbox.com/s/p0hlgk5jot1ipzc/stackd.cache.3.png?raw=1)

This chart shows the relative frequency (in %) of each possible outcome on the network shown above, for the {0,0,1} cache configuration (i.e. a request has to go up to tier 3 to get to the content). The x axis shows all possible FP and ALPHA combinations. 

By the way, the possible outcomes are:

* **Correct:** Request is delivered to the correct cache, i.e. C<sub>3</sub>.
* **Incorrect:** Request is delivered to an incorrect destination (due to the occurrence FPs). In this case, this never happens.
* **Feedback/Fallback:** Request is relayed to the origin server (which in this case coincides with C<sub>3</sub>). There are two strategies for relays:
	* **Feedback:** Say that a FP is detected at router in t2. The request is sent back to the source, which then sends it towards C<sub>3</sub>.
	* **Fallback:** If a FP is detected at a router, it is immediately sent towards C<sub>3</sub> and not sent back to the requester. This results in lower latencies (see ).
* **Dropped:** If a router doesn't know what to do with a request, it simply drops the packet.

### Avg. Latencies
<a name="subsubsec:avg-lat"></a>

![](https://www.dropbox.com/s/auxh8j6p11fnela/bar.cache.2.png?raw=1)

This one shows how avg. latencies vary for each FP & ALPHA combination. This is for the case in which both C<sub>2</sub> and C<sub>3</sub> can correctly serve the request. Therefore, the minimum latency is 5 hops (green line): 2 hops to get up to t2, 1 hop to peer to another domain in t2, 3 hops down to C<sub>2</sub>. The origin server is considered to be at C<sub>3</sub>: by the same reasoning, its latency can be easily derived to be 7 hops (red line). Obviously, the 'feedback' (pink) relaying method takes more time than 'fallback' (light green).

**How do we calculate this avg. value?** Each outcome can happen with a certain probability. Our model calculates the probability for each possible outcome, given a particular <FP per tier, ALPHA per tier, cache distance> configuration (editable via `.scn` files in `test/configs/example/`. Also, each one of these outcomes will have an associated latency (in nr. of hops). All we have to do in the end is:

$\text{avg. latency}=\sum^{\forall O_i} P_i\,L_i$

, in which:

* $O_i$: Some outcome $i$, e.g. 'request packet delivered to a particular incorrect destination through t2'.
* $P_i$: Probability of $O_i$ happening. $\sum^{\forall i} P_i = 1$.
* $L_i$: Latency for outcome $O_i$ (in nr. of hops).

### Probability DAGs
<a name="subsubsec:prob-dag"></a>

If you want, you can visualize what the possible outcomes are, along with their associated probabilities, in the form of a probability Directed Acyclic Graph (DAG). We use a tool called [Graphviz](http://www.graphviz.org/) for that. 

E.g. to do so for the case of high FP rates (`hfp`) and low ALPHA (`la`), you have to run:

`$ cd test/data/example/cache.2/`

`$ dot -Tpdf hfp.la.dot -o ../../../graphs/example/hfp.la.cache.2.pdf`

`$ open ../../../graphs/example/hfp.la.cache.2.pdf`

This will yield the following DAG:

![](https://www.dropbox.com/s/d33ck2jsn7jumfr/hfp.la.cache.2.png?raw=1)

Here's how you roughly read it:

* `ORI00`: Root of the DAG, the initial outcome (i.e. 'the request is issued by R<sub>1</sub>').
* `<outcome>(<curr. level>:<to level>)`: A possible router lookup outcome (intermediate). An RID router can experience the following lookup outcomes when looking up an RID in its forwarding engine:
	* `TPO`: Lookup **O**nly returns TRUE postive entries. The packet is correctly forwarded.
	* `SFP`: Lookup yields a **S**ingle match, and it's a FALSE positive. The packet is incorrectly forwarded.
	* `MHS`: Lookup yields multiple matches (or **H**its), all pointing to the **S**ame interface. This could include both TRUE and FALSE positives. If there are TRUE positives in the table, the packet is correctly forwarded; if not, the packet is incorrectly forwarded.
	* `MHD`: Lookup yields multiple matches, pointing to the **D** interface. A FP is detected mid-way. The packet is immediately relayed.
	* `DEF`: Lookup cannot find a positive entry in the RID table, either FALSE or TRUE. Packet is either sent to the tier above OR dropped (if we're already at the top tier). 
* `<outcome>(<curr. level>:EOP)`: Similar to the above, but in this case, the packet reaches reaches the end of its life as an RID packet ('**E**nd **O**f **P**ath'). This can happen if the packet is delivered to a correct/incorrect destination, or if it is relayed upon the occurrence of an `MHD` lookup outcome at some router.
* Each edge has the associated probability of jumping between 2 lookup particular outcomes.