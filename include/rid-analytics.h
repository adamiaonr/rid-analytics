#ifndef RID_ANALYTICS_HH
#define RID_ANALYTICS_HH

#define MAX_FILENAME_SIZE   64
#define MAX_REQUEST_SIZE    30
#define MAX_TIER_DEPTH      30

#define END_OF_PATH         (int) -1

// interface events
#define EVENT_NUM       0x04 // hack : we don't count w/ EVENT_TTL
#define EVENT_NLM       0x00 // no link matches
#define EVENT_MLM       0x01 // multiple link matches
#define EVENT_LLM       0x02 // local link match
#define EVENT_SLM       0x03 // single link match (other than local)
#define EVENT_TTL       0x04 // drop due to rtt expiration
#define EVENT_UNKNOWN   0x05

#define RES_FALLBACKS       0x04

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <sstream>
#include <fstream>  // reading .nw files

#include "tree/tree.hh"
#include "rid-router.h"
#include "path-state.h"

typedef std::map<std::string, std::vector<uint8_t> > RID_TPMap;
typedef std::map<std::string, RID_Router *> RID_RouterMap;
typedef std::pair<std::string, RID_Router *> RID_RouterPair;

class RID_Analytics {

    public:

        struct run_record
        {
            RID_Router * next_router;
            int in_iface;
            Path_State * state;
        };

        RID_Analytics() {}
        RID_Analytics(    
            std::string nw_filename,
            uint8_t request_size,
            uint16_t bf_size,
            std::string origin_server_id,
            int mode,
            std::string output_dir,
            std::string output_label);
        ~RID_Analytics();

        void run(std::string scn_filename, std::string start_router_id);

    private:

        int build_network(std::string scn_filename);
        int read_scn(std::string scn_filename);
        int run_rec(
                RID_Router * router, 
                uint8_t ingress_iface,
                Path_State * prev_state);
        int erase_access_tree_rec(RID_Router * router);

        int get_origin_distance(RID_Router * from_router);
        int get_origin_distance_rec(RID_Router * from_router);

        int handle_nlm(RID_Router * router, __float080 event_prob, Path_State * prev_state);
        int handle_llm(RID_Router * router, __float080 event_prob, Path_State * prev_state);
        int handle_mlm(RID_Router * router, __float080 event_prob, Path_State * prev_state);
        int handle_slm(
            RID_Router * router,
            Path_State * prev_state,
            uint8_t in_iface,
            std::vector<std::vector<__float080> > iface_probs,
            std::vector<std::vector<__float080> > out_fptree_probs,
            std::vector<RID_Analytics::run_record> & slm_stack);
        int handle_ttl_drop(
            RID_Router * router,
            __float080 event_num,
            Path_State * prev_state);

        void init_output(std::string label, std::string output_dir);
        void add_outcome(RID_Router * router, int outcome, __float080 prob, __float080 latency);
        void add_events(RID_Router * router, std::vector<__float080> event_probs);
        void add_fwd(RID_Router * router, int in_iface, Path_State * state, std::vector<std::vector<__float080> > iface_probs);

        uint8_t request_size;
        uint16_t bf_size;
        int mode;

        std::vector<__float080> f_r_distribution;

        // for quick access to RID_Router * in the topology, via id strings 
        RID_TPMap tp_sizes;
        RID_RouterMap routers;
        // AS associated w/ origin server
        RID_Router * origin_server;
        RID_Router * start_router;

        // ttl for network topology
        int ttl;

        // array of output files (w/ results)
        std::vector<std::ofstream> output_file;
};

#endif