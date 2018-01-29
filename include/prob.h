#ifndef PROB_HH
#define PROB_HH

// FIXME: tried to use __float128, didn't work... instead of writing 'long 
// double', decided to go with the alias '__float080' (80 bit size)
#define __float080 long double

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>

class Prob {

    public:

    struct fp_data {

        uint8_t tp_size = 0;                // max tp size (note: this initialization is permitted by C++11)
        __float080 num_entries = 0.0;       // nr. of entries
        __float080 entry_prop = 0.0;        // proportion of entries
        std::vector<__float080> f_distr;    // distr. of entry sizes
        std::vector<__float080> f_r_distr;  // |f\r| distr.

        bool on_fptree = false;             // is iface on a fp tree? FIXME: should have size info
        bool is_blocked = false;            // is the iface blocked?
    };

    Prob(
        __float080 m,
        __float080 n,
        uint8_t iface_num) {

        // bf size and req. size (as bf fp rate params)
        this->m = m;
        this->n = n;
        // k : optimal nr. of hash functions for bfs
        this->k = (log(2) * ((__float080) this->m)) / ((__float080) this->n);

        this->iface_num = iface_num;
        // for (uint8_t f = 0; f < this->n; f++)
        //     this->iface_on_fptree_prob.push_back(0.0);

        // initialize largest match prob matrix
        for (unsigned i = 0; i < this->iface_num; i++) {
            this->lm_prob.push_back(std::vector<std::vector<__float080> > ((n + 1), std::vector<__float080> ((n + 1), 0.0)));
            this->anti_lm_prob.push_back(std::vector<std::vector<__float080> > ((n + 1), std::vector<__float080> ((n + 1), 0.0)));
        }
        // mark largest match probs as 'not calc'ed'
        this->has_lm_prob = false;

        // initialize P('iface i on fp tree of size p') matrix
        this->iface_on_fptree_prob = std::vector<std::vector<__float080> >((iface_num), std::vector<__float080> ((n + 1), 0.0));
    }

    ~Prob() {}

    int calc_iface_probs(
        std::vector<std::vector<fp_data> > * iface_fp_data,
        std::vector<__float080> * in_fptree_prob,
        std::vector<__float080> & iface_probs,
        bool strict = false);

    private:

    int calc_lm_probs(std::vector<std::vector<fp_data> > * iface_fp_data);
    int calc_lm_prob(fp_data iface_fp_data, uint8_t i, bool anti = false);
    int calc_iface_on_fptree_probs(
        std::vector<std::vector<fp_data> > * iface_fp_data,
        std::vector<__float080> * in_fptree_prob);

    __float080 calc_iface_prob(uint8_t i, bool strict = false);

    void calc_log_prob_fp_neq(Prob::fp_data iface_fp_data, std::vector<__float080> & log_prob_fp_neq);
    void print_lm_prob(uint8_t iface, bool anti = false);

    // largest match probabilities
    std::vector<std::vector<std::vector<__float080> > > lm_prob;
    std::vector<std::vector<std::vector<__float080> > > anti_lm_prob;
    // P('iface i on fp tree of size p') matrix
    std::vector<std::vector<__float080> > iface_on_fptree_prob;

    // fp rate parameters
    __float080 m;   // fp rate m : nr. of bits in bf
    __float080 n;   // fp rate n : nr. of elements encoded in bf (same as req size)
    __float080 k;   // fp rate k : nr. of hashes in bf
    // num of ifaces in respective router
    // FIXME: why do we need this?
    uint8_t iface_num;

    // flag to mark if lm probs are cached or not
    bool has_lm_prob;
};

#endif