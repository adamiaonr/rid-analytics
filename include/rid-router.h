#ifndef RID_ROUTER_HH
#define RID_ROUTER_HH

// FIXME: tried to use __float128, didn't work... instead of writing 'long 
// double', decided to go with the alias '__float080' (80 bit size)
#define __float080 long double

// various max. limits
#define MAX_RID_ROUTER_STRING_SIZE  128
#define MAX_RID_ROUTER_IFACE_NUM    10
#define MAX_RID_ROUTER_ENTRY_SIZE   30

// special iface indexes
#define IFACE_LOCAL     0x00

// interface events
#define EVENT_NUM       0x04 // hack : we don't count w/ EVENT_RTT
#define EVENT_NLM       0x00 // no link matches
#define EVENT_MLM       0x01 // multiple link matches
#define EVENT_LLM       0x02 // local link match
#define EVENT_SLM       0x03 // single link match (other than local)
#define EVENT_RTT       0x04 // drop due to rtt expiration
#define EVENT_UNKNOWN   0x05

// modes for calc_cumulative_prob()
#define MODE_EI_EXCLUSIVE   0x00
#define MODE_MIS            0x01
#define MODE_LI             0x02
#define MODE_EI_INCLUSIVE   0x03

// modes for handling multiple matches in routers. this 
// influences efficiency. there are 3 diff. modes:
//  -# MMH_FLOOD    : forward over all matching ifaces
//  -# MMH_RANDOM   : forward over 1 random iface
//  -# MMH_FALLBACK : forward using the fallback address
#define MMH_FLOOD           0x00
#define MMH_RANDOM          0x01
#define MMH_FALLBACK        0x02

// modes for egress prefix tree probability calculation
#define EGRESS_PTREE_PROB_MODE_LOCAL    0x00
#define EGRESS_PTREE_PROB_MODE_GLOBAL   0x01

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

class RID_Router {

    public:

        struct nw_address
        {
            RID_Router * router;
            uint8_t iface;
        };

        struct lpm_pmf_row
        {
            __float080 * lpm_pmf_prob;
        };

        struct fwd_table_row
        {
            // interface nr. (0 is always the local iface)
            uint8_t iface;
            // % of entries in fwd_table pointing to iface
            __float080 iface_proportion;        
            // distribution of fwd entry sizes of iface. e.g. for IFACE_LOCAL 
            // one could expect a distr. skewed towards f = f_max, while other 
            // ifaces a distr. skewed towards f = 1 (single prefixes like 
            // '/cmu')
            __float080 * f_distribution;
            // bitmask of fwd entry sizes. if iface has > 0 entries w/ size 
            // f, the f-th bit of f_bitmask is set to 1, otherwise it's set 
            // to 0.
            uint16_t f_bitmask;
            // nr. of prefixes associated with this iface
            uint64_t num_entries;
            // bitmasks of source trees included in this iface
            uint8_t * tree_bitmask;
            int tree_bitmask_size;
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
            uint8_t f_max,
            uint16_t bf_size,
            int mm_mode);
        ~RID_Router();

        int init(
            std::string router_id,
            uint64_t fwd_table_size,
            uint8_t iface_num,
            uint8_t f_max,
            uint16_t bf_size,
            int mm_mode);

        int add_fwd_table_entry(
            int iface, 
            __float080 iface_proportion, 
            std::map<int, __float080> * size_dist,
            std::vector<uint8_t> * tree_bitmask,
            RID_Router * next_hop_router,
            int next_hop_iface);

        void set_fwd_table(RID_Router::fwd_table_row * fwd_table);
        RID_Router::fwd_table_row * get_fwd_table();
        uint64_t get_fwd_table_size();
        int get_num_entries(uint8_t iface);
        RID_Router::nw_address get_next_hop(uint8_t iface);

        std::string get_id() { return this->id; }
        uint8_t get_iface_num();
        uint8_t get_f_max();

        int forward(
            uint8_t request_size,
            uint8_t ingress_iface,
            int * tp_sizes,
            uint8_t * tree_bitmask,
            __float080 ingress_prob,
            __float080 * ingress_ptree_prob,
            __float080 * f_r_distribution);

        void print_lpm_pmf(uint8_t iface);

        __float080 get_iface_events_prob(uint8_t event);
        void print_iface_events_pmf();

        __float080 get_egress_iface_prob(uint8_t iface);
        __float080 * get_egress_iface_probs(uint8_t iface);
        void print_egress_iface_prob();

        int get_tree_bitmask_size(uint8_t iface);
        uint8_t * get_tree_bitmask(uint8_t iface);

        std::set<uint8_t> get_blocked_ifaces() {
            return this->no_forwarding;
        }

    private:

        // computation of LPM match distributions (L_{i,p})
        void init_lpm_pmf(uint8_t iface);
        __float080 * get_lpm_prob(uint8_t iface, uint8_t ptree_size);
        __float080 get_lpm_prob(uint8_t iface, uint8_t ptree_size, uint8_t f);
        int calc_largest_match_distributions(
            uint8_t request_size, 
            int * tp_sizes, 
            __float080 * f_r_distribution);

        // computation of joint distribution of L_{i,p}, for all ifaces
        // void init_joint_lpm_pmf(__float080 ** joint_prob_matrix);
        // void clear_joint_lpm_pmf(__float080 ** joint_prob_matrix);
        void add_joint_lpm_prob(std::map<std::string, __float080> * joint_prob_matrix, int * iface_pivots, __float080 value);
        __float080 get_joint_lpm_prob(std::map<std::string, __float080> * joint_prob_matrix, int * iface_pivots);
        int calc_joint_largest_match_distributions(
            uint8_t ptree_size, 
            uint8_t ptree_iface,
            std::map<std::string, __float080> * joint_prob_matrix);

        // computation of iface event & egress size probabilities
        int calc_iface_events_distributions(
            int ptree_size,
            std::map<std::string, __float080> * joint_prob_matrix);

        int calc_log_prob_not_fp(
            __float080 m, 
            __float080 n, 
            __float080 k,
            __float080 f_entries,
            __float080 * f_distribution, 
            __float080 * f_r_distribution,
            __float080 * log_prob_not_fp);

        __float080 calc_cumulative_prob(
            uint8_t mode,
            uint8_t fixed_iface, 
            uint8_t fixed_iface_size,
            std::map<std::string, __float080> * joint_prob_matrix);

        __float080 calc_log_joint_largest_match_prob(
            uint8_t ptree_size,
            uint8_t ptree_iface, 
            int * iface_pivots);

        __float080 calc_joint_largest_match_prob_sum(
            uint8_t ptree_size,
            uint8_t ptree_iface);

        bool is_iface_on_ptree(int iface, uint8_t * tree_bitmask);
        int calc_ptree_iface_probs();

        RID_Router::lpm_pmf_row ** lpm_pmf;
        // probability of interface events (NIS, MIS, LI)
        __float080 * iface_events_pmf;
        // probability of having the largest match at iface i
        __float080 ** egress_iface_prob;
        // ingress probabilities, what are these?
        //  -# ingress_prob : the probability of having the request forwarded 
        //                    INTO this router (by the previous router), regardless 
        //                    of the cause. the causes can be 3:
        //                      -# true positive (in which case ingress_prob = 1.0)
        //                      -# stuck in a 'false positive' prefix tree. the 
        //                         prefixes announced in the tree can be of 
        //                         size 1, 2, 3, ..., f_max. the probability for a 
        //                         specific size is given by the ingress_ptree_prob 
        //                         array (see below).
        //                      -# fallback address
        //
        //  -# ingress_ptree_prob : the probability of having the request 'stuck' in
        //                          a 'false positive' prefix tree of size ptree_size,
        //                          with 1 <= ptree_size <= f_max.
        __float080 ingress_prob;
        __float080 * ingress_ptree_prob;
        // router id follows the format <tree_index>.<height>.<width>
        std::string id;
        // number of interfaces (at least 2 : LOCAL & // UPSTREAM)    
        uint8_t iface_num;
        // keep track of interfaces over which a router CAN'T forward packets 
        // (e.g. ingress iface, upstream ifaces for valley-free routing, etc.)
        std::set<uint8_t> no_forwarding;
        // max. possible size of a forwarding entry 
        uint8_t f_max;
        uint16_t bf_size;              
        // size of forwarding table (important for some calculations)
        uint64_t fwd_table_size;
        // the forwarding table : an array of fwd_table_row structs
        RID_Router::fwd_table_row * fwd_table;
        // if the RID router is initialized
        bool initialized;
        // multiple match handling mode
        uint8_t mm_mode;

        // tells if iface i is potentially a continuation of a prefix tree 
        std::vector<bool> iface_on_ptree;
        std::vector<__float080> ptree_iface_probs;
};

#endif