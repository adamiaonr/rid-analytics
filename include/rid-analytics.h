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

// indexes for output file array
#define FILE_EVENTS         0
#define FILE_PATHS          1

// error handling mode
#define EH_DEFAULT      0x00
#define EH_FALLBACK     0x01

// P(|F|_i) distributions (for IFACE_LOCAL and otherwise)
#define LOCAL       0x00
#define NON_LOCAL   0x01

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "tree/tree.hh"
#include "rid-router.h"
#include "path-state.h"

typedef std::map<std::string, int *> RID_TPMap;
typedef std::map<std::string, RID_Router *> RID_RouterMap;
typedef std::pair<std::string, RID_Router *> RID_RouterPair;

class RID_Analytics {

    public:

        struct nw_link
        {
            int iface_local;
            int iface_remote;
            std::string router_remote;
        };

        RID_Analytics() {}
        RID_Analytics(
            std::string nw_filename,
            uint8_t request_size,
            uint16_t bf_size,
            std::string origin_server,
            int mm_mode,
            int eh_mode);
        ~RID_Analytics();

        int run(std::string scn_filename);
        int view_results(uint8_t mode, std::string output_dir, std::string output_label);

    private:

        int build_network(std::string scn_filename);
        int read_scn(std::string scn_filename);
        int run_rec(
                RID_Router * router, 
                uint8_t ingress_iface,
                tree<Path_State *>::iterator prev_path_state_itr);
        int erase_access_tree_rec(RID_Router * router);

        int get_origin_distance(RID_Router * from_router);
        int get_origin_distance_rec(RID_Router * from_router);
        int on_path_to_origin(RID_Router * router, int iface);

        // NETWORK PARAMETERS : 

        // height of tree
        uint8_t access_tree_height;
        // max. possible size for requests & forwarding entries
        uint8_t f_max;

        // FORWARDING PARAMETERS : only necessary when calling forward() on an 
        // RID router
        uint8_t request_size;
        uint16_t bf_size;
        int mm_mode;
        int eh_mode;
        RID_TPMap tp_sizes;
        __float080 * f_r_distribution;
        // path state tree
        tree<Path_State *> path_state_tree;
        // for quick access to RID_Router * in the topology, via id strings 
        RID_RouterMap routers;
        // AS associated w/ origin server
        RID_Router * origin_server;
        // AS path to origin server
        std::set<std::string> origin_server_path;
};

#endif