<a name="sec:intro"></a>
# rid-analytics

Implementation of an analytical model of a forwarding method based on [Bloom filters](https://en.wikipedia.org/wiki/Bloom_filter). This forwarding method is designated as 'RID forwarding' in the [eXpressive Internet Architecture](http://www.cs.cmu.edu/~xia/resources/Documents/XIA-nsdi.pdf) (XIA) (check [this repo](https://github.com/adamiaonr/xia-core/tree/xia-v2-rids) for an implementation), and was adapted from an initial proposal published in [Papalini et al. 2014](http://nis-ita.org/ITA_static/attachments/2791/icn8945.pdf).

`rid-analytics`' uderlying engine is based on a probabilistic graphical model which determines the probability of having an RID packet follow any given path in a network topology (e.g. PoP-level topologies made available by the [Rocketfuel ISP mapping engine](http://research.cs.washington.edu/networking/rocketfuel/)). Details about the math behind it are available in [this paper](https://www.dropbox.com/s/cjixlvjrbhtbjl2/infocom-2018-extended.pdf?dl=0).

<a name="sec:usage"></a>
# Compilation

```bash
$ git clone https://github.com/adamiaonr/rid-analytics.git
$ cd rid-analytics
$ make
```
<a name="sec:usage"></a>
# Usage

There are 2 ways to use `rid-analytics`:
1. Direct use of `run-analysis` binary (not reccommended)
2. Using `.test` files

<a name="subsec:direct-run"></a>
## Direct use of `run-analysis` binary

The compilation produces the basic binary of `rid-analytics`, which runs individual experiments: `run-analysis`. It can be used as follows:

```bash
$ ./run-analysis --scn-file <path to .scn file> --data-dir <path to results directory> --output-label <label to prepend to result files> --bf-size <bloom filter size in bit> --request-size <nr. of url elements in request> --mm-mode <multiple match res. code> --resolution-mode <enable/disable delivery error resolution> --origin-server <id of origin server> --start-router <id of router which issues request>
```

It requires multiple inputs, of which the most important is an `.scn` file. `.scn` files describe a experiment to run, including: 
* Router adjacencies in the topology under evaluation
* Distribution of # of entries per {router, link} pair
* Distribution of entry sizes per {router, link} pair
* Shortest paths from router *i* every other router *j*

Fortunately, `.scn` files can be generated automatically, using the `generate_tests.py` script in `scripts/rocketfuel/`. Another important input is `--data-dir`, on which the `.tsv` files with the results of the experiments will be saved. Result files are explained in (...).

<a name="subsec:test-files"></a>
## Using `.test` files

`rid-analytics` comes with a set of Python scripts which automatically create experiment scenarios, based on PoP-level topologies from [Rocketfuel](http://research.cs.washington.edu/networking/rocketfuel/). The flow is as follows:
1. Use `scripts/rocketfuel/generate_tests.py` to generate a `.test` file and a batch of `.scn` files (referred by the `.test` file)
2. Use `scripts/rocketfuel/run_evaluation.py` to run the batch of experiments described in the `.test` file

<a name="subsubsec:generate-test-files"></a>
### Generating `.test` and `.scn` files

As an example, we generate a batch of experiments with the following specs:
* **Topology:** 4755 VSNL (India)
* **Request sizes:** 5, 10
* **Bloom filter sizes:** 192, 256
* **Entry sizes:** all entries should be size 1
* **Fwd. table size:** 10M entries
* **Multiple match resolution:** Random Matched Link (RML) forwarding
* **Delivery error resolution:** Feedback
* **# of source/destination pairs:** 3 pairs, 4 hops away from each other
* **True positives:** 1 true positive, size 1, 4 hops away from requester

To do that, run:
```bash
$ cd <work dir>/rid-analytics/scripts/rocketfuel
$ python generate_tests.py --test-dir <work dir>/rid-analytics/experiments/examples/rocketfuel/tests --topology-file <work dir>/rid-analytics/experiments/examples/rocketfuel/pop-level-maps/4755/edges.wt --req-sizes 5:10 --bf-sizes 192:256 --entry-sizes "1:100" --table-sizes 10000000 --modes 0:0 --path-sizes 3:4 --add-tps 1:1:4
```
After this, the directory specified in the `--test-dir` option should contain 3 folders: 
1. `configs`: filled w/ `.test` and `.scn` files 
2. `results`: empty 
3. `topologies`: filled w/ `.pdf` files depicting the chosen topologies, in graph form (example of [4755 VSNL (India) topology](https://www.dropbox.com/s/lg99ab6h4ogzl8u/infocom-2018.pdf?dl=0))

<a name="subsubsec:running-experiments"></a>
### Running experiments

Run:
```bash
$ cd <work dir>/rid-analytics/scripts/rocketfuel
$ python run_test.py --test-file <work dir>/rid-analytics/experiments/examples/rocketfuel/tests/configs/4755.test
```
This should fill the `<work dir>/rid-analytics/experiments/examples/rocketfuel/tests/` directory with a bunch of `.tsv` files, with the results from the experiments.

<a name="subsec:result-files"></a>
## Result files

Result files are saved in `.tsv` format. There are 2 types of files: **event** files and **path** files.

<a name="subsubsec:event-files"></a>
### Event files

Event files report the probability of witnessing specific **forwarding events** at some router, for all routers which can be visited by a request packet during an experiment. There are 5 types of forwarding events:
0. No link matches (NLM)
1. Multiple link matches (MLM)
2. Local link matches (LLM)
3. Single link match (SLM)
4. TTL drop (TTL)

***FIXME:** Understanding the probabilities of these events is a bit confusing, since these are not mutually exclusive.*

Following the running example, the first lines of the file `<work dir>/rid-analytics/experiments/examples/rocketfuel/tests/results/events.4755-192-05-01-10000000-0-0-6-4.1501589180.tsv` read:

```
ROUTER  EVENT	PROB
6	0	0
6	1	0
6	2	0
6	3	1
3	0	0
3	1	0.0818887
3	2	0.02097
3	3	0.918111
0	0	0
0	1	0.0414116
0	2	0.0310071
0	3	0.958588

(...)
```

Let's use a graph depiction of the topology to make this easier to follow...

![](https://www.dropbox.com/s/r452hf3ezxlk74b/4755.png?raw=1)

Starting at router 6, there is a 100% likelihood of having a single link match, on the link to router 3. On router 3, there's a ~0.02 chance of having a cache hit (LLM event), and a ~0.92 chance of having a single link match towards router 0. There's a ~0.08 chance of having multiple link matches: in this case the joint event LLM AND a match on the link towards router 0.

<a name="subsubsec:path-files"></a>
### Path files

Path files report the evolution of **path states** as a request transitions among routers during an experiment. Similarly to event files, path files report the probability of having some state **s** at each router in the path. There are 5 types of path states:
0. Correct delivery
1. Incorrect delivery
2. Fallback delivery
3. Fallback relay
4. Packet drop
5. TTL drop

Following the running example, the first lines of the file `<work dir>/rid-analytics/experiments/examples/rocketfuel/tests/results/path.4755-192-05-01-10000000-0-0-6-4.1501589180.tsv` read:

```
ROUTER	STATUS	LATENCY	PROB
6	4	0	0
6	6	0	0
6	1	0	0
3	4	1	0
3	6	1	0
3	1	1	0.02097
0	4	2	0
0	6	2	0
0	1	2	0.0310071
1	4	3	0
1	6	3	0
1	1	3	0.0107375

(...)
```

Again, starting at router 6: as there's no way of terminating a path on router 6 (remember that event SLM, i.e. of having the packet pass to the next router, has probability 1.0), there's no probability associated with any state at that point. At the next router - router 3 - there's some probability of having the path end in an *incorrect delivery*, with probability equal to ~0.02 (equal to the probability of event LLM, local link match). The latency at that point is 1 hop, corresponding to the hop in-between routers 6 and 3.
