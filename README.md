<a name="sec:intro"></a>
# rid-analytics

Implementation of an analytical model of a forwarding method based on [Bloom filters](https://en.wikipedia.org/wiki/Bloom_filter). This forwarding method is designated as 'RID forwarding' in the [eXpressive Internet Architecture](http://www.cs.cmu.edu/~xia/resources/Documents/XIA-nsdi.pdf) (XIA) (check [this repo](https://github.com/adamiaonr/xia-core/tree/xia-v2-rids) for an implementation), and was adapted from an initial proposal published in [Papalini et al. 2014](http://nis-ita.org/ITA_static/attachments/2791/icn8945.pdf).

`rid-analytics`' uderlying engine is based on a probabilistic graphical model which determines the probability of having an RID packet follow any given path in a network topology (e.g. PoP-level topologies made available by the [Rocketfuel ISP mapping engine](http://research.cs.washington.edu/networking/rocketfuel/)).

<a name="sec:usage"></a>
# Usage

<a name="subsec:scn"></a>
## Example

<a name="subsec:run"></a>
## Scenario file generation

## Running the analytics engine on scenario files

## Output files and visualization

<a name="sec:analytical-model"></a>
# Analytical Model

<a name="subsec:analytical-model-overview"></a>
## Overview

At its core, the analytical model of RID forwarding calculates the probability of specific forwarding decisions within a router, assuming Longest Prefix Matching (LPM) semantics. E.g. say a router has 3 interfaces - 0, 1 and 2 - each with associated forwarding entries which can go from sizes 1 to *p*. At its core, our model calculates the probability of events such as 'having the router decide to send the packet over interface 1 due to a match of size 3'.

To better illustrate how this works, we introduce a simple example, depicted in the diagram below. 

The main characteristics of this example are:
* A topology with 5 routers: R0 to R4
* Each router holds forwarding entries of size 1 or 2 
* Imagine each router is the gateway to a lower-level access network which advertises 1 million prefixes (represented as caches in the diagram). As such, each router has 1 million forwarding entries associated with the local interface (interface 0), and 4 million entries distributed over the remaining interfaces.
* Following from above, the forwarding table at each router holds a total of 5 million entries
* The flow of the packet - say a request for some piece of content - starts at router R0
* The requested content is held by router R2 **only**. Other routers don't hold the requested content. Let's assume R2 announces that content using a prefix of size 1.
* The shortest paths from R0 to any other router are represented by green dotted lines

<a name="subsec:analytical-model-step-1"></a>
## Step 1 : Forwarding at the initial router (R0)

### Forwarding events

We're interested in calculating the probabilities of forwarding events at router R0. We can represent such events as an array with size 3 - one index for each interface - with each position capable of holding values from 0 to 2, indicating the size of the *largest* match at the interface associated with that index. E.g. the array \[0,1,1\] represents the event of having no matches at interface 0, and the largest match of size 1 at interfaces 1 and 2.

Let's list the events which can happen at the initial router. These events will ultimately decide the fate of a packet flow in the topology:
1. Since the shortest path from R0 to R2 goes through interface 1, the probability of having a match of size 1 over interface 1 is 1.0 (i.e. certain). In our model, this constraint is modeled by the presence of a true positive (TP) entry of size 1 at interface 1.
2. Following from the above, it is impossible to have no match at any of the interfaces (since a match of size 1 is guaranteed at interface 1). In other words, the probability of event \[0,0,0\] is impossible: in fact, any event like \[\*,0,\*\] is impossible, since index 1 must at least have value 1. 
3. A false positive match of size *p* (with *p* equal to 1 or 2) happens at interface *i*. Note that this includes having more than 1 interface experiencing a match of size *p* simultaneously.
4. Note that the occurrence of FP matches at an interface also determine the probability of having a packet `bind' to a false positive tree of size *p*. E.g. if the event is \[0,1,1\] happens at the first router, the packet will surely bind to a FP tree of size 1 starting at interface 2. In addition, if a FP match of size 1 happened at interface 1 (together with the TP match), the packet will also bind to a FP tree of size 1 at interface 1.

To calculate the probability of each of the different events listed above, we follow the algorithm below. We will explain each step in detail in the subsequent sub-sections.

```
// probability of having a match of size m on iface i, when iface i is part of 
// a tree of size p
largest_match_prob[len(router.get_ifaces())][max_fwd_entry_size][max_tree_size]

// hash table of joint event probabilities. event keys follow a specific 
// string format. e.g. for event [0,1,1], the key is '000101'. only events 
// with probability > 0.0 are added to the table.
joint_largest_match_prob[event_key_str]

// probability of having the packet forwarded over iface i due to a match of 
// size m
iface_probs[len(router.get_ifaces())][max_fwd_entry_size]

for iface in router.get_ifaces():
    
    // if iface doesn't have any associated entries,
    // add it to a 'no forwarding' list and skip any 
    // iface processing
    if iface.get_number_of_entries() == 0:
        invalid_ifaces.add(iface)
        continue

    // calculate the probabilities of having matches 
    // of any size (0, 1 or 2) for iface. this fills 
    // largest_match_prob[iface].
    calc_largest_match_probs(iface)

// calculate the event probabilities for the case of having 
// a packet bound to a prefix tree of size tree_size
for tree_size in [0, ..., max_fwd_entry_size]:

    // if the probability of having a packet enter the router 
    // while bound to a tree of size tree_size is 0.0, then 
    // the contribution of any related events to final event 
    // probabilities is also 0.0. as such, skip any subsequent calculations. 
    if prob_ingress_tree[tree_size] == 0.0:
        continue

    // calculate the contribution to the final event probabilities 
    // for the case of having the prefix tree of size ptree_size on 
    // iface
    for iface in router.get_ifaces():

        if iface in invalid_ifaces:
            continue

        // calculate the probability of all events, for the 
        // case of having the prefix tree of size ptree_size 
        // continue over iface, and add this contribution to 
        // joint_largest_match_prob
        add_joint_largest_match_probs(tree_size, iface)

// finally, by scanning joint_largest_match_prob, calculate the 
// cumulative probabilities which give the probability of having 
// the router forward a packet over iface, with a match of size m
calc_iface_probs()

return iface_probs
```

### Largest match distributions per interface

In this step, we calculate the probability of having a packet match a forwarding entry of size *m*, at some interface *i*. We represent the match size per interface as a random variable *L<sub>i</sub>*, In our example, this step should output an array of size 3 for each of the interfaces, one value for each of the possible match lengths: 0 (i.e. no match, 1 and 2).

Calculation of $P(L_{i} = m)$ is based on the basic formulas for [Bloom filter false positive rates](https://en.wikipedia.org/wiki/Bloom_filter), namely

![equation](http://www.sciweavers.org/tex2img.php?eq=1%2Bsin%28mc%5E2%29&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=)


### *P(L<sub>i,p</sub>, ..., L<sub>j,p</sub>)*: Probability of forwarding events

### *P(L<sub>i,p</sub>, ..., L<sub>j,p</sub>)*: Probability of FP tree bindings

<a name="subsec:analytical-model-step-2"></a>
## Step 2 : Forwarding at router R1

After router R0, the packet arrives at router R1 (it can also arrive at R4, but we skip that step in this description). There, we calculate the same probabilities listed above. Since R1 is not the initial router, there are some differences in the calculations, due to the effect FP tree bindings. We explain them below.

### *P(L<sub>i,p</sub>, ..., L<sub>j,p</sub>)*: Probability of forwarding events w/ FP tree bindings

