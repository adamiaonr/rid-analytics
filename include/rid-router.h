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

// interface events
#define EVENT_NUM           0x04

#define EVENT_NIS           0x00
#define EVENT_MIS           0x01
#define EVENT_LI            0x02
#define EVENT_EI            0x03
//#define EVENT_UNKNOWN       0x04

// modes for calc_cumulative_prob()
#define MODE_EI             0x00
#define MODE_MIS            0x01
#define MODE_LI             0x02

// default values for Bloom Filters
#define DEFAULT_M   (__float080) 160.0

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>

class RID_Router {

    public:

        struct lpm_pmf_row
        {
            // // L_{i,p} = {0, 1, ..., |R|_max} for iface i
            // uint8_t iface;

            // // ifaces can be associated with a 'prefix tree', to which the 
            // // request got 'stuck'. this is the size of the 
            // // tree (\in {1, ..., |R|_max}).
            // uint8_t ptree_size;
            // 
            __float080 * lpm_pmf_prob;
            // __float080 * lpm_not_in_ptree_pmf;
        };

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

        void set_as_starting_router();
        bool is_starting_router();

        int forward(
            uint8_t request_size,
            uint8_t ingress_iface,                   
            int * tp_sizes,     
            __float080 * ingress_probs,          
            __float080 * f_r_distribution);

        void print_lpm_pmf(uint8_t iface);

        __float080 get_iface_events_prob(uint8_t event);
        void print_iface_events_pmf();

        __float080 get_egress_size_prob(uint8_t iface, uint8_t f);
        void print_egress_size_pmf();

    private:

        // computation of LPM match distributions (L_{i,p})
        void init_lpm_pmf(uint8_t iface);
        __float080 * get_lpm_prob(uint8_t iface, uint8_t ptree_size);
        __float080 get_lpm_prob(uint8_t iface, uint8_t ptree_size, uint8_t f);
        int calc_lpm_pmf(
            uint8_t request_size, 
            uint8_t ingress_iface, 
            int * tp_sizes, 
            __float080 * f_r_distribution);

        // computation of joint distribution of L_{i,p}, for all ifaces
        void init_joint_lpm_pmf(__float080 ** joint_prob_matrix);
        void clear_joint_lpm_pmf(__float080 ** joint_prob_matrix);
        void add_joint_lpm_prob(__float080 * joint_prob_matrix, int * iface_pivots, __float080 value);
        __float080 get_joint_lpm_prob(__float080 * joint_prob_matrix, int * iface_pivots);
        int calc_joint_lpm_pmf(__float080 * joint_prob_matrix, uint8_t ptree_size, uint8_t ptree_iface);

        // computation of iface event & egress size probabilities
        void init_egress_size_pmf(uint8_t iface);
        int calc_iface_events_pmf(__float080 * joint_prob_matrix);

        int get_log_fp_rates(
            __float080 m, 
            __float080 n, 
            __float080 k,
            __float080 iface_proportion,
            __float080 * f_distribution, 
            __float080 * f_r_distribution,
            __float080 * log_fp_rates);

        __float080 calc_cumulative_prob(__float080 * joint_prob_matrix, uint8_t iface, uint8_t f, uint8_t mode);
        __float080 calc_joint_log_prob(uint8_t ptree_size, uint8_t in_ptree_iface, int * iface_pivots);

        // pmf for random variables L_i : size of the longest 
        // match pointing to iface i. L_i can assume values in 
        // {0, 1, 2, ..., f_max}.
        RID_Router::lpm_pmf_row ** lpm_pmf;

        // // joint pmf of random variables L_i (L_0 x L_1 x ... x L_n)
        // // FIXME: note that if L_0 > 0, path automatically ends at the 
        // // router
        // __float080 * joint_lpm_size_pmf;

        // pmf for random variable I, the interface to be chosen by the LPM 
        // engine. can take values {-2, -1, 0, 1, 2, ..., iface_num - 1}:
        //  * -2 is used for multiple iface hits
        //  * -1 is used for no hits (remember 0 is the IFACE_LOCAL)

        __float080 total_joint_prob;

        // probability of interface events (NIS, MIS, LI)
        __float080 * iface_events_pmf;

        // ingress size probabilities
        __float080 * ingress_size_pmf;
        // probability of having iface i output associated with the longest 
        // match of size L
        __float080 ** egress_size_pmf;

        // a scenario is composed by n access trees
        uint8_t access_tree_index;  
        // each router as (x, y) coordinates in the access tree
        uint8_t height;             
        uint8_t width;
        // if this is the starting router, set it to TRUE
        bool starting_router;

        // number of interfaces (at least 2 : LOCAL & // UPSTREAM)    
        uint8_t iface_num;          
        uint8_t iface_ingress;
        // max. possible size of a forwarding entry 
        uint8_t f_max;              
        // size of forwarding table (important for some calculations)
        uint32_t fwd_table_size;

        // the forwarding table : an array of fwd_table_row structs
        RID_Router::fwd_table_row * fwd_table;
};

#endif