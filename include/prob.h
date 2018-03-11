#ifndef PROB_HH
#define PROB_HH

// FIXME: tried to use __float128, didn't work... instead of writing 'long 
// double', decided to go with the alias '__float080' (80 bit size)
#define __float080 long double

// interface events
#define EVENT_NUM       0x04 // hack : we don't count w/ EVENT_TTL
#define EVENT_NLM       0x00 // no link matches
#define EVENT_MLM       0x01 // multiple link matches
#define EVENT_LLM       0x02 // local link match
#define EVENT_SLM       0x03 // single link match (other than local)
#define EVENT_TTL       0x04 // drop due to rtt expiration
#define EVENT_UNKNOWN   0x05

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

    class fp_data {

        public:

        fp_data() {

            tp_size = 0;
            num_entries = 0.0;
            entry_prop = 0.0;

            on_fptree = false;
            is_blocked = false;
        }
        ~fp_data() {}

        // FIXME: to make things easy, let them be public
        uint8_t tp_size;                    // max tp size (note: this initialization is permitted by C++11)
        __float080 num_entries;             // nr. of entries
        __float080 entry_prop;              // proportion of entries
        std::vector<__float080> f_distr;    // distr. of entry sizes
        std::vector<__float080> f_r_distr;  // |f\r| distr.

        std::vector<uint8_t> tree_bitmask;  // bitmask of iface
        bool on_fptree;
        bool is_blocked;                    // is the iface blocked?
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

        // initialize prob arrays:
        //  - P(L_i = l | P_i = p)
        //  - P(L_(~i) = l | P_(~i) = p)
        for (unsigned i = 0; i < this->iface_num; i++) {
            this->lm_cond_prob.push_back(std::vector<std::vector<__float080> > ((n + 1), std::vector<__float080> ((n + 1), 0.0)));
            this->lm_complement_cond_prob.push_back(std::vector<std::vector<__float080> > ((n + 1), std::vector<__float080> ((n + 1), 0.0)));
        }
        //  - P(L_i = l)
        lm_marg_prob = std::vector<std::vector<__float080> > ((this->iface_num), std::vector<__float080>((this->n + 1), 0.0));
        //  - P(L_{~i} = l)
        lm_complement_marg_prob = std::vector<std::vector<__float080> > ((this->iface_num), std::vector<__float080>((this->n + 1), 0.0));
        //  - P(P_i = p) and P(P_(~i) = p)
        this->fptree_prob.push_back(std::vector<std::vector<__float080> > ((this->iface_num), std::vector<__float080> ((n + 1), 0.0)));
        this->fptree_prob.push_back(std::vector<std::vector<__float080> > ((this->iface_num), std::vector<__float080> ((n + 1), 0.0)));

        // mark largest match probs as 'not calc'ed'
        this->has_lm_prob = false;
    }

    ~Prob() {}

    int calc_probs(
        std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
        std::vector<uint8_t> * tree_bitmask,    // tree bitmask from prev router
        __float080 in_prob,
        std::vector<__float080> * in_fptree_prob,
        std::vector<std::vector<__float080> > & iface_probs,
        std::vector<__float080> & event_num,
        std::vector<std::vector<__float080> > & out_fptree_probs);

    private:

    int calc_lm_probs(
        std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
        std::vector<uint8_t> * tree_bitmask,
        std::vector<std::vector<__float080> > & out_fptree_probs);
    int calc_lm_prob(
        uint8_t i, 
        Prob::fp_data * iface_fp_data,
        std::vector<uint8_t> * tree_bitmask,
        std::vector<std::vector<__float080> > & out_fptree_probs, 
        bool iface_complement = false);
    int calc_iface_prob(
        uint8_t i, 
        std::vector<std::vector<Prob::fp_data *> > * iface_fp_data, 
        std::vector<__float080> & iface_prob);

    int calc_event_num(
        std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
        std::vector<std::vector<__float080> > iface_probs,
        std::vector<__float080> & event_num);

    int calc_fptree_probs(
        std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
        std::vector<uint8_t> * tree_bitmask,
        std::vector<__float080> * in_fptree_prob);

    int calc_out_fptree_prob(
        uint8_t i,
        Prob::fp_data * iface_fp_data,
        std::vector<__float080> log_prob_fp_neq,
        std::vector<__float080> log_prob_fp_smeq,
        std::vector<std::vector<__float080> > & out_fptree_probs);

    void calc_log_prob_fp_neq(
        Prob::fp_data * iface_fp_data,
        std::vector<__float080> & log_prob_fp_neq,
        __float080 lt_ratio = 1.0);
    void print_lm_prob(uint8_t i, bool iface_complement = false);

    // basic units of prob calculation:
    //  - P(L_i = l | P_i = p)
    std::vector<std::vector<std::vector<__float080> > > lm_cond_prob;
    //  - P(L_(~i) = l | P_(~i) = p)
    std::vector<std::vector<std::vector<__float080> > > lm_complement_cond_prob;
    //  - P(L_i = l)
    std::vector<std::vector<__float080> > lm_marg_prob;
    //  - P(L_{~i} = l)
    std::vector<std::vector<__float080> > lm_complement_marg_prob;
    //  - P(P_i = p) and P(P_(~i) = p)
    std::vector<std::vector<std::vector<__float080> > > fptree_prob;

    // fp rate parameters
    __float080 m;   // fp rate m : nr. of bits in bf
    __float080 n;   // fp rate n : nr. of elements encoded in bf (same as req size)
    __float080 k;   // fp rate k : nr. of hashes in bf

    // num of ifaces in respective router
    // FIXME: why do we need this?
    uint8_t iface_num;
    // flag to mark if lm probs are cached or not
    // FIXME: why do we need this?
    bool has_lm_prob;
};

#endif