<a name="sec:intro"></a>
# rid-analytics

Implementation of an analytical model of a forwarding method based on [Bloom filters](https://en.wikipedia.org/wiki/Bloom_filter). This forwarding method is designated as 'RID forwarding' in the [eXpressive Internet Architecture](http://www.cs.cmu.edu/~xia/resources/Documents/XIA-nsdi.pdf) (XIA) (check [this repo](https://github.com/adamiaonr/xia-core/tree/xia-v2-rids) for an implementation), and was adapted from an initial proposal published in [Papalini et al. 2014](http://nis-ita.org/ITA_static/attachments/2791/icn8945.pdf).

`rid-analytics`' uderlying engine is based on a probabilistic graphical model which determines the probability of having an RID packet follow any given path in a network topology (e.g. PoP-level topologies made available by the [Rocketfuel ISP mapping engine](http://research.cs.washington.edu/networking/rocketfuel/)). Details about the math behind it are available in [this paper](https://www.dropbox.com/s/cjixlvjrbhtbjl2/infocom-2018-extended.pdf?dl=1).

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

## Direct use of `run-analysis` binary

The compilation produces the basic binary of `rid-analytics`, which runs individual experiments: `run-analysis`. It can be used as follows:

```bash
$ ./run-analysis --scn-file <path to .scn file> --data-dir <path to results directory> --output-label <label to prepend to result files> --bf-size <bloom filter size in bit> --request-size <nr. of url elements in request> --mm-mode <multiple match res. code> --resolution-mode <enable/disable delivery error resolution> --origin-server <id of origin server> --start-router <id of router which issues request>
```

It requires multiple inputs, of which the most important is an `.scn` file. `.scn` files describe a experiment to run, including: 
* Router adjacencies in the topology under evaluation
* Distribution of # of entries per {router, link}
* Distribution of entry sizes per {router, link}
* Shortest paths from router *i* every other router *j*

Fortunately, `.scn` files can be generated autmatically, using the `generate_tests.py` tool in `scripts/rocketfuel`. Another important input is `--data-dir`, on which the `.tsv` files with the results of the experiments will be saved. Result files are explained in (...).

## Using `.test` files

Instead of 
