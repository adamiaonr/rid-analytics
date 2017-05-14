<a name="sec:intro"></a>
# rid-analytics

Implementation of an analytical model of a forwarding method based on [Bloom filters](https://en.wikipedia.org/wiki/Bloom_filter). This forwarding method is designated as 'RID forwarding' in the [eXpressive Internet Architecture](http://www.cs.cmu.edu/~xia/resources/Documents/XIA-nsdi.pdf) (XIA) (check [this repo](https://github.com/adamiaonr/xia-core/tree/xia-v2-rids) for an implementation), and was adapted from an initial proposal published in [Papalini et al. 2014](http://nis-ita.org/ITA_static/attachments/2791/icn8945.pdf).

`rid-analytics`' uderlying engine is based on a probabilistic graphical model which determines the probability of having an RID packet follow any given path in a network topology (e.g. PoP-level topologies made available by the [Rocketfuel ISP mapping engine](http://research.cs.washington.edu/networking/rocketfuel/)).

<a name="sec:analytical-model"></a>
# Analytical Model

At its core, the analytical model of RID forwarding calculates the probability of specific forwarding decisions within a router, assuming Longest Prefix Matching (LPM) semantics. E.g. say a router has 3 interfaces - 0, 1 and 2 - each with associated forwarding entries with sizes 1 or 2. The model calculates the probability of events such as 'having a packet forwarded over interface 1 due to a match of size 1'.

To better illustrate how this works, we introduce a running example, depicted in the diagram below.

![](https://www.dropbox.com/s/2c9y0882dcejtma/model-example-crop.png?raw=1)

The main characteristics of this example are:
* A topology with 5 routers: R0 to R4
* Each router holds forwarding entries of size 1 or 2 
* Imagine each router is the gateway to a lower-level access network which advertises 1 million prefixes (represented as caches in the diagram). As such, each router has 1 million forwarding entries associated with the local interface (interface 0), and 4 million entries distributed over the remaining interfaces.
* Following from above, the forwarding table at each router holds a total of 5 million entries
* The flow of the packet - say a request for some piece of content - starts at router R0
* The requested content is held by router R2 **only**. Other routers don't hold the requested content. Let's assume R2 announces that content using a prefix of size 1.
* The shortest paths from R0 to any other router are represented by green dotted lines

## Main algorithm

To calculate the probability of forwarding decisions, we follow the algorithm below. We will explain each step in detail in the next sub-sections.

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
    
    if iface.get_number_of_entries() == 0:
        invalid_ifaces.add(iface)
        continue

    calc_largest_match_probs(iface)

// calculate the event probabilities for the case of having 
// a packet bound to a prefix tree of size tree_size
for tree_size in [0, ..., max_fwd_entry_size]:

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

calc_iface_probs()

return iface_probs
```

## Forwarding events

### Objective and notation

Calculate the probability of *forwarding events* at a router. Regarding router R0 in our example, we represent such events as an 1 x 3 array - one index per interface - each value indicating the size of the *largest* match at the interface associated with that index. E.g. the array \[0,1,1\] represents the event of having no matches at interface 0, and a largest match of size 1 at interfaces 1 and 2. The probability of this forwarding event is *P(*\[0,1,1\]*)*.

### Causes of forwarding events

Forwarding events have 4 underlying causes:
1. **True positive (TP) matches:** Since the shortest path from R0 to R2 goes through interface 1, the probability of having a match of size 1 over interface 1 is 1.0 (i.e. certain). In our model, this constraint is modeled by the presence of a TP entry of size 1 at interface 1.
2. **True negative (TN) matches:** Following from the above, it is impossible to have no matches at any of the interfaces (since a match of size 1 is guaranteed at interface 1). In other words, the probability of event \[0,0,0\] is impossible: in fact, any event like \[\*,0,\*\] is impossible, since index 1 must at least have value 1. 
3. **False positive (FP) matches:** A FP match of size *p* (with *p* equal to 1 or 2) happens at interface *i*. Note that this includes having more than 1 interface experiencing a match of size *p* simultaneously.
4. **Matches due to FP prefix tree bindings:** At any point in its path, a packet may bind to a FP prefix tree of size *p*. The reasoning is: if a FP match of size 1 happens at R0, then the same FP match will happen at subsequent routers, since the forwarding entry which causes the match at R0 is announced on a tree sourced at a subsequent router. If a packet is bound to a prefix tree of size *p*, then a match of size *p* is guaranteed to happen. If the packet is not bound to a tree, the occurrence of a 'new' FP match at an interface determines the probability of having the packet `bind' to a new FP prefix tree of size *p*. E.g. if the event is \[0,1,1\] happens at R0, the packet will bind to a FP tree of size 1 starting at interface 2. Likewise, if a FP match of size 1 happened at interface 1 (together with the TP match), the packet will also bind to a FP tree of size 1 at interface 1.

The diagram below shows the possible prefix trees which a packet starting at R0 can bind itself to.

![](https://www.dropbox.com/s/jlh87l2iu4t64i6/model-example-ptrees-crop.png?raw=1)

## Largest match probabilities per interface

### Objective and notation

Calculate the probability of having *m* as the largest positive match size on interface *i*, assuming a packet bound to a prefix tree of size *p*. We represent the match size per interface as a random variable *L<sub>i,p</sub>*, In our example, this step should output a 3 x 3 array per interface, as shown below:

<table>
  <tr>
    <td> </td>
    <td colspan="3">Largest match size <i>m</i></td>
  </tr>
  <tr>
    <td>Prefix tree size <i>p</i></td>
    <td>0</td>
    <td>1</td>
    <td>2</td>
  </tr>
  <tr>
    <td>0</td>
    <td><i>P(L<sub>i,0</sub> = 0)</i></td>
    <td><i>P(L<sub>i,0</sub> = 1)</i></td>
    <td><i>P(L<sub>i,0</sub> = 2)</i></td>
  </tr>
  <tr>
    <td>1</td>
    <td><i>P(L<sub>i,1</sub> = 0)</i></td>
    <td>...</td>
    <td>...</i></td>
  </tr>
  <tr>
    <td>2</td>
    <td><i>P(L<sub>i,2</sub> = 0)</i></td>
    <td>...</td>
    <td>...</i></td>
  </tr>
</table>

### Calculation details

The calculation of *P(L<sub>i,p</sub>)* is as follows:

![](https://www.dropbox.com/s/5iw4k7kdptwlh2o/eq-calc-largest-probs.jpg?raw=1)

The reasoning behind the branches of the above expression is (top to bottom):
1. If interface *i* has a true positive (TP) entry larger than *m* (i.e. |*TP*| *> m*), then the largest match will always be *at least* |*TP*|. As such, the probability of having the largest match equal to *m* on interface *i* is 0. The same applies if *i* is the continuation of a prefix tree of size larger than *m*.
2. Even if interface *i* has a TP entry of size *m*, *m =* |*TP*| will only be the largest match size iif there is **not** a FP match larger than |*TP*|.
3. For match sizes larger than the TP entry size *m >* |*TP*| size, we calculate the intersection of 2 events (we assume the independence of the events): having a FP match of size *m* **and** not having a FP match larger than *m*.

| Symbol                              | Description                                |
| ----------------------------------- |:------------------------------------------:|
| m                                   | Bloom filter size (in bit) |
| \|*R*\| or \|*F*\|<sub>*max*</sub>  | Request size or max. forwarding entry size |
| *T*                                 | Table size in nr. of entries               |
| *T(i)*                              | Nr. of entries per interface               |
| *P(*\|*F*\|,*i)*                    | Distribution of entry sizes per interface  |
| *P(*\|*F*\\*R*\|,*i)*               | Distribution of request-entry differences per interface |

**Basic FP rate expressions:** The FP probabilities - e.g. *P(*|*FP*|*= m)* - are based on the basic expressions for [Bloom filter false positive rates](https://en.wikipedia.org/wiki/Bloom_filter), given below, and dependent on the parameters listed in the Table above.


![](https://www.dropbox.com/s/1bd1evislmrnwwk/total-router-level.jpg?raw=1)

![](https://www.dropbox.com/s/9ly2xjed4hhqftc/simple-fpr.jpg?raw=1)


## Joint probability of largest matches (forwarding events)

### Objective and notation

Calculate the probability of joint largest match events, i.e. the probability of having a largest match of size m at interface i, size n at interface j, and so on. This is essentially the probability of forwarding events such as \[0,1,1\] in R0 in our example. The calculation of *P(*\[*m<sub>0</sub>*,*m<sub>1</sub>*,...,*m<sub>i</sub>*\]*)* is a bit complex, involving 2 steps: (1) the intersection of largest match probabilities for each interface; and (2) the contribution of FP prefix trees.

### Probability calculation

**Intersection of largest match probabilities:** Assuming that largest match events are independent in-between interfaces, the intersection of largest match probabilities - e.g. *P(*\[0,1,1\]*)* - can be calculated as *P(*\[0,1,1\]*) = P(L<sub>0,p</sub> = 0) * P(L<sub>1,p</sub> = 1) * P(L<sub>2,p</sub> = 1)*. The complex part is in handling the contribution of FP prefix trees.

**Contributions of FP prefix trees:** Say the packet arrives from router R0 into R1's interface 1, with 0.6 probability of being bound to a FP prefix tree of size 1. This can happen when a FP match of size 1 happens in parallel with the expected TP match of size 1 at R0's interface 1. We hereby refer to such tree as an *ingress prefix tree* of size 1, and represent its probability as *P(*|*ptree*|*<sub>ingress</sub> = 1) = 0.6* (in the section [Probability of egress prefix trees](#probability-of-egress-prefix-trees), we show how the ingress prefix tree probabilities (*P(*|*ptree*|*<sub>ingress</sub>)*) are calculated).

Then, with 0.6 probability, that ingress tree must continue on one of R1's outgoing interfaces (i.e. 0, 2 **or** 3). This implies that a match of size 1 is guaranteed at that interface. Note that we did not include interface 4 as a possible continuation for an ingress tree. This is because the link R1:R4 is not in the list of shortest paths announced into R0. As such, any FP match which results into a FP tree binding at R0's interface 1 must come from either R2 or R3. In any case, a new FP prefix tree can start at interface 4.

**General expression:** More generally, say we want to calculate the probability of event \[*m<sub>0</sub>*,*m<sub>1</sub>*,*m<sub>2</sub>*,*m<sub>3</sub>*\] at R1. Regarding the association of prefix trees to interfaces, this event can happen if (1) the prefix tree continues on interface 0; (2) if it continues on interface 2; **or** (3) if it continues on interface 3. We assume these events are independent, i.e. if the prefix tree continues over interface 0, then it won't continue over interface 2 or 3. As such, the probability of the forwarding event is:

![](https://www.dropbox.com/s/i1f6ga5hj7m74tj/forwarding-event-prob.jpg?raw=1)

The second term within the final is the calculation of the joint largest match probability in-between interfaces. We first consider the probability of a largest match of size *m<sub>i</sub>* on interface *i*, assuming *i* to be the continuation of the prefix tree of size *p*. Then, the largest match probabilities for the remaining interfaces assume those interfaces to not be associated with any prefix tree (hence *p = 0*). Finally, to assure that the sum of forwarding event probabilities is 1.0, we normalize the second term.

**Probability of interfaces and prefix trees:** In general, the probability of having interface *i* as the continuation of the ingress tree of size *p* - *P(*|*ptree*|*<sub>i</sub> = p)* - is based on the fraction of the announcements of size *p* between the egress interface of the previous router - R0's interface 0 in our example - and a possible egress interface *i*:

![](https://www.dropbox.com/s/1u32vzh3b6p19lq/real_ptree_prob.jpg?raw=1)

To determine the numerator and denominator of the above expression, we keep the information about reachable sources in `.scn` files (source bitmasks per interface). 

However, note that at the initial router - R0 - there is no 'previous router', since the packet starts its path there. As such, at R0, *P(*|*ptree*|*<sub>ingress</sub> = 0) = 1* and *P(*|*ptree*|*<sub>i</sub> = p) = (1 / (# of egress ifaces))*, i.e. the binding probabilities are equally distributed equally over the egress interfaces.

## Probability of interface choices

### Objective and notation

Calculate the probability of having packet forwarded over egress interface *i*, due to a match of size *m*. We use the notation *P(I = i, M = m)*. At router R0 in our example, interface choice events are represented as a 3 x 3 array (3 interfaces, 3 possible match sizes), as shown below:

<table>
  <tr>
    <td> </td>
    <td colspan="3">Largest match size <i>m</i></td>
  </tr>
  <tr>
    <td>Interface <i>i</i></td>
    <td>0</td>
    <td>1</td>
    <td>2</td>
  </tr>
  <tr>
    <td>0</td>
    <td><i>P(I = 0, M = 0)</i></td>
    <td><i>P(I = 0, M = 1)</i></td>
    <td><i>P(I = 0, M = 2)</i></td>
  </tr>
  <tr>
    <td>1</td>
    <td><i>P(I = 1, M = 0)</i></td>
    <td><i>...</i></td>
    <td><i>...</i></td>
  </tr>
  <tr>
    <td>2</td>
    <td><i>P(I = 2, M = 0)</i></td>
    <td><i>...</i></td>
    <td><i>...</i></td>
  </tr>
</table>

### Probability calculation

*P(I = i, M = m)* can be calculated using the forwarding event probabilities, in 2 modes:
1. **Exclusive mode:** Probability of having the largest match of size *m* at interface *i* only, i.e. other interfaces should have largest matches **smaller** than *m*. The probability becomes *P(I = i, M = m) = P(*\[*m<sub>0</sub> < m*,*m<sub>1</sub> < m*,...,*m<sub>i</sub> = m*\]*)*
2. **Inclusive mode:** Probability of having the largest match of size *m* at interface *i*, regardless of also having other largest matches of size *m* at other interfaces. The probability becomes *P(I = i, M = m) = P(*\[*m<sub>0</sub> <= m*,*m<sub>1</sub> <= m*,...,*m<sub>i</sub> = m*\]*)*

## Probability of egress prefix trees

### Objective and notation

Finally, calculate the probability of FP tree bindings for each of the egress interfaces. These will then serve as ingress prefix tree probabilities (*P(*|*ptree*|*<sub>ingress</sub>) = m*) for the routers connected to each of the egress interfaces. We use the notation *P(*|*ptree*|*<sub>egress i</sub>) = m* to refer to the egress prefix tree probabilities of interface *i*. At R0 in our example, we represent these probabilities as a 3 x 3 array (3 interfaces, 3 possible prefix tree sizes).

### Probability calculation

An interface *i* can be the continuation of a prefix tree of size *m* iif:
1. There are common source announcements in-between interface *i* and the announcements on the egress interface of a previous router. [In our example](https://www.dropbox.com/s/jlh87l2iu4t64i6/model-example-ptrees-crop.png?raw=1), R1's interface 3 has a common announcement with those coming in R1's interface 1 ('R3' AND \['R1', 'R2', 'R3'\] == 'R3'). 
2. There are no common announcements, but the probability of a FP match of size *m* is larger than 0.0
3. There is a TP entry on interface *i*, and as such, even if there is no other match besides a TP match, the packet continues to the next router. In this case, interface *i* is determined to be in the prefix tree of size 0.

The cases above translate into the following expression:

![](https://www.dropbox.com/s/3qe2b1sdm9zp8z0/egress-ptree-probs.jpg?raw=1)

As a final step, we normalize the egress prefix tree probabilities as follows:

![](https://www.dropbox.com/s/34iq1nugq0ar61k/egress-ptree-probs-norm.jpg?raw=1)
