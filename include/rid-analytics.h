#ifndef RID_ANALYTICS_HH
#define RID_ANALYTICS_HH

#define MAX_FILENAME_SIZE   64
#define MAX_REQUEST_SIZE    30
#define MAX_TIER_DEPTH      30

#define END_OF_PATH         (int) -1

#define AVERAGE_LATENCY     (char *) "average_latency"
#define CACHE_LATENCY       (char *) "cache_latency"
#define ORIGIN_LATENCY      (char *) "origin_latency"

#define PENALTY_TYPE_FEEDBACK   (int) 0
#define PENALTY_TYPE_FALLBACK   (int) 1

#define PENALTY_TYPE_STR_SIZE   2

#define MODE_VERBOSE        0x01
#define MODE_SAVE_CDF       0x02
#define MODE_SAVE_GRAPH     0x04
#define MODE_SAVE_OUTCOMES  0x08

// P(|F|_i) distributions (for IFACE_LOCAL and otherwise)
#define LOCAL       0x00
#define NON_LOCAL   0x01

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "tree/tree.hh"
#include "rid-router.h"
#include "path-state.h"

class RID_Analytics {

    public:

        RID_Analytics() {}
        RID_Analytics(
            uint8_t access_tree_num,
            uint8_t access_tree_height,
            uint8_t iface_num,
            uint8_t f_max,
            int f_min_annc,
            int f_min,
            int expand_factor,
            uint16_t bf_size,
            uint64_t fwd_table_size,
            __float080 * iface_entry_proportion,
            __float080 ** f_distributions);
        ~RID_Analytics();

        int print_tp_sizes();
        int print_iface_entry_proportions();

        int run(
            uint8_t request_size,
            int * tp_sizes,
            __float080 * f_r_distribution);

        int view_results(uint8_t mode, char * output_file_path);

    private:

        int build_network();
        int build_network_rec(RID_Router * parent_router);
        int run_rec(
                RID_Router * router, 
                uint8_t ingress_iface,
                tree<Path_State *>::iterator prev_path_state_itr);
        int erase_access_tree_rec(RID_Router * router);

        int get_tp_distance(RID_Router * from_router);
        int get_tp_distance_rec(RID_Router * from_router);

        // NETWORK PARAMETERS : 

        // # of of access trees in the network
        uint8_t access_tree_num;
        // height and outdegree of the tree (outdegree is iface_num - 1)
        uint8_t access_tree_height;
        uint8_t iface_num;
        // max. possible size for requests & forwarding entries
        uint8_t f_max;
        int f_min_annc;
        // min. possible size for forwarding entries
        int f_min;
        int expand_factor;
        // % of entries in pointing to an iface
        __float080 * iface_entry_proportion;
        // distribution of |F| over iface (i.e. % of entries 
        // with |F| = x). e.g. for IFACE_LOCAL one could expect a 
        // distr. skewed towards |F| = f_max, while for IFACE_UPSTREAM a 
        // distr. skewed towards |F| = 1 (single prefixes like '/cmu')
        __float080 * f_distribution_local;
        __float080 * f_distribution_non_local;
        // forwarding table size
        uint64_t fwd_table_size;

        // an array of access trees (pointers to root of access trees)
        RID_Router ** root_routers;
        // special pointer to 'starting' router of the network (by convention, 
        // bottom left router of access tree at index 0)
        RID_Router * start_router;

        // FORWARDING PARAMETERS : only necessary when calling forward() on an 
        // RID router
        uint8_t request_size;
        uint16_t bf_size;
        int * tp_sizes;
        __float080 * f_r_distribution;

        // path state tree
        tree<Path_State *> path_state_tree;
};

#endif