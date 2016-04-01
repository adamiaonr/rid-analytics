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
#define IFACE_MULTI_HITS    -2
#define IFACE_NO_HITS       -1
#define IFACE_LOCAL         0
#define IFACE_UPSTREAM      1
#define IFACE_DOWNSTREAM    2

// default values for Bloom Filters
#define DEFAULT_M   (__float080) 160.0

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "path-state.h"

class RID_Router {

    public:

        struct fwd_table_row
        {
            // interface (LOCAL, UPSTREAM, DOWNSTREAM, etc.)
            // FIXME: this info might be redundant
            int iface;      

            // % of entries in fwd_table pointing to iface
            __float080 iface_proportion;        
            // distribution of |F| over iface (i.e. proportion of entries 
            // with |F| = x). e.g. for IFACE_LOCAL one could expect a 
            // distr. skewed towards |F| = f_max, while for IFACE_UPSTREAM a 
            // distr. skewed towards |F| = 1 (single prefixes like '/cmu')
            __float080 * f_distribution;

            // the RID_Router object at the other end of this iface
            RID_Router * next_hop;
        };

        RID_Router() {}
        RID_Router(
            uint8_t access_tree_index, 
            uint8_t height, 
            uint8_t width,
            uint32_t fwd_table_size,
            uint8_t iface_num,
            uint8_t f_max);
        ~RID_Router();

        int add_fwd_table_entry(
            int iface, 
            __float080 iface_proportion, 
            __float080 * f_distribution);
        void set_fwd_table_next_hop(uint8_t iface, RID_Router * next_hop_router);
        RID_Router * get_fwd_table_next_hop(uint8_t iface);
        void set_fwd_table(RID_Router::fwd_table_row * fwd_table);
        RID_Router::fwd_table_row * get_fwd_table();
        uint32_t get_fwd_table_size();

        uint8_t get_iface_num();
        uint8_t get_f_max();
        uint8_t get_access_tree_index();
        uint8_t get_height();
        uint8_t get_width();

        int forward(
            uint8_t request_size,
            uint8_t ingress_iface,                   
            int * tp_sizes,                 
            __float080 * f_r_distribution);

        void print_f_pmf();
        void print_joint_f_pmf();

        __float080 get_lpm_iface_pmf(int iface);
        void print_lpm_iface_pmf();

    protected:

        __float080 * init_f_pmf(uint8_t iface);
        void set_f_pmf(int iface, int f, __float080 prob);
        __float080 * get_f_pmf(uint8_t iface);
        __float080 get_f_pmf(uint8_t iface, uint8_t f);
        int calc_f_pmf(
            uint8_t request_size,
            uint8_t ingress_iface,
            int * tp_sizes, 
            __float080 * f_r_distribution);

        __float080 * init_joint_f_pmf();
        void set_joint_f_pmf(int * iface_pivots, __float080 value);
        __float080 get_joint_f_pmf(int * iface_pivots);
        int calc_joint_f_pmf();

        __float080 * init_lpm_iface_pmf();
        void set_lpm_iface_pmf(int iface, __float080 value);
        int calc_lpm_iface_pmf();

    private:

        int get_log_fp_rates(
            __float080 m, 
            __float080 n, 
            __float080 k,
            __float080 iface_proportion,
            __float080 * f_distribution, 
            __float080 * f_r_distribution,
            __float080 * log_fp_rates);

        __float080 calc_cumulative_joint_f(uint8_t iface, uint8_t f);
        __float080 calc_joint_f_log_prob(int * iface_pivots);
        __float080 calc_no_match_prob();

        // pmf for the random variable F_i, the size of the longest match 
        // pointing to iface i. F_i can assume values in {0, 1, 2, ..., f_max}.
        __float080 ** f_pmf;

        // joint pmf of random variables F_i (F_0, F_1, ..., F_n)
        __float080 * joint_f_pmf;

        // pmf for random variable I, the interface to be chosen by the LPM 
        // engine. can take values {-2, -1, 0, 1, 2, ..., iface_num - 1}:
        //  * -2 is used for multiple iface hits
        //  * -1 is used for no hits (remember 0 is the IFACE_LOCAL)
        __float080 * lpm_iface_pmf;

        // a scenario is composed by n access trees
        uint8_t access_tree_index;  
        // each router as (x, y) coordinates in the access tree
        uint8_t height;             
        uint8_t width;

        // size of forwarding table (important for some calculations)
        uint32_t fwd_table_size;
        // number of interfaces (at least 2 : LOCAL & // UPSTREAM)    
        uint8_t iface_num;          
        // max. possible size of a forwarding entry 
        uint8_t f_max;              

        // the forwarding table : an array of fwd_table_row structs
        RID_Router::fwd_table_row * fwd_table;
};

#endif