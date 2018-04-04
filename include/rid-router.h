#ifndef RID_ROUTER_HH
#define RID_ROUTER_HH

// FIXME: tried to use __float128, didn't work... instead of writing 'long 
// double', decided to go with the alias '__float080' (80 bit size)
#define __float080 long double

// various max. limits
#define MAX_RID_ROUTER_STRING_SIZE  128
#define MAX_RID_ROUTER_IFACE_NUM    10
#define MAX_RID_ROUTER_ENTRY_SIZE   30

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include <string.h>

#include <map>
#include <set>
#include <algorithm>
#include <iostream>
#include <sstream>      // std::ostringstream

#include "prob.h"

class RID_Router {

    public:

        struct nw_address {

            RID_Router * router;
            uint8_t iface;
        };

        struct fwd_table_row {

            // interface nr. (0 is always the local iface)
            uint8_t iface = -1;

            // % of entries in fwd_table pointing to iface
            __float080 iface_proportion;        
            // distribution of fwd entry sizes of iface. e.g. for IFACE_LOCAL 
            // one could expect a distr. skewed towards f = f_max, while other 
            // ifaces a distr. skewed towards f = 1 (single prefixes like 
            // '/cmu')
            std::vector<__float080> f_distr;
            std::vector<__float080> f_r_distr;
            // bitmask of fwd entry sizes. if iface has > 0 entries w/ size 
            // f, the f-th bit of f_bitmask is set to 1, otherwise it's set 
            // to 0.
            uint16_t f_bitmask;
            // nr. of prefixes associated with this iface
            uint64_t num_entries;
            // bitmasks of source trees included in this iface
            std::map<int, std::vector<uint8_t> > tree_bitmasks;
            // pointer to an RID router (including the ingress iface on the 
            // other end)
            RID_Router::nw_address next_hop;
        };

        RID_Router() {
            this->initialized = false;
        };

        RID_Router(
            std::string router_id,
            uint64_t fwd_table_size,
            uint8_t iface_num,
            uint8_t req_size,
            uint16_t bf_size,
            int mode);

        ~RID_Router();

        int init(
            std::string router_id,
            uint64_t fwd_table_size,
            uint8_t iface_num,
            uint8_t req_size,
            uint16_t bf_size,
            int mode);

        std::string get_id() { return this->id; }

        uint8_t get_req_size() { return this->req_size; }
        uint8_t get_bf_size() { return this->bf_size; }
        uint8_t get_iface_num() { return this->iface_num; }
        uint64_t get_fwd_table_size() { return this->fwd_table_size; }
        uint64_t get_num_entries(uint8_t i) {return (int) this->fwd_table[i].num_entries; }

        // add forwarding entry info
        int add_fwd_table_entry(
            int iface, 
            __float080 iface_proportion, 
            std::map<int, __float080> * f_distr,
            std::map<int, std::vector<uint8_t> > * tree_bitmask,
            RID_Router * next_hop_router,
            int next_hop_iface);

        RID_Router::nw_address get_next_hop(uint8_t i) { return this->fwd_table[i].next_hop; }
                
        int get_tree_bitmask_size(uint8_t i) { return this->fwd_table[i].tree_bitmasks[0].size(); }
        std::vector<uint8_t> * get_tree_bitmask(uint8_t i) { return &(this->fwd_table[i].tree_bitmasks[0]); }
        std::map<int, std::vector<uint8_t> > * get_tree_bitmasks(uint8_t i) { return &(this->fwd_table[i].tree_bitmasks); };

        std::vector<bool> * get_blocked_ifaces() { return &(this->blocked_ifaces); }

        std::vector<uint8_t> * get_tp_sizes() { return &(this->tp_sizes); };
        void set_tp_sizes(std::vector<uint8_t> * tp_sizes) { this->tp_sizes = (*tp_sizes); };
        void set_f_r_distr(std::vector<__float080> * f_r_distr) {
            for (uint8_t i = 0; i < this->iface_num; i++)
                this->fwd_table[i].f_r_distr = (*f_r_distr);
        };

        int forward(
            uint8_t ingress_iface,
            std::map<int, std::vector<uint8_t> > * tree_bitmasks,
            __float080 ingress_prob,
            std::vector<__float080> * in_fptree_prob,
            std::vector<std::vector<__float080> > & iface_probs,
            std::vector<std::vector<__float080> > & out_fptree_probs);

    private:

        void is_iface_on_fptree(int i, std::map<int, std::vector<uint8_t> > * tree_bitmasks);
        Prob::fp_data * get_fp_data(RID_Router * router, uint8_t i, bool iface_complement = false);

        // router id follows the format <tree_index>.<height>.<width>
        std::string id;
        // if the RID router is initialized
        bool initialized;
        // bloom filter parameters
        uint8_t req_size;
        uint16_t bf_size;
        // number of ifaces (at least 2 : LOCAL & // UPSTREAM)    
        uint8_t iface_num;
        // size of forwarding table (important for some calculations)
        uint64_t fwd_table_size;
        // the forwarding table : an array of fwd_table_row structs
        std::vector<RID_Router::fwd_table_row> fwd_table;
        // tp sizes for this router
        std::vector<uint8_t> tp_sizes;
        // keep track of interfaces over which a router CAN'T forward packets
        std::vector<bool> blocked_ifaces;
        // is iface i a continuation of a prefix tree of size t?
        std::vector< std::vector<bool> > iface_on_fptree;
        // probability module
        Prob * prob_mod;
};

#endif