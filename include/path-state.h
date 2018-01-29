#ifndef PATH_STATE_HH
#define PATH_STATE_HH

#define MAX_PATH_STATE_STRING_SIZE 128

// possible outcomes
#define OUTCOME_CORRECT_DELIVERY    (int) 0x00
#define OUTCOME_INCORRECT_DELIVERY  (int) 0x01
#define OUTCOME_FALLBACK_DELIVERY   (int) 0x02
#define OUTCOME_FALLBACK_RELAY      (int) 0x03
#define OUTCOME_PACKET_DROP         (int) 0x04
#define OUTCOME_TTL_DROP            (int) 0x05
#define OUTCOME_UNDEF               (int) 0x06

// path status
#define STATUS_TP     (int) 0x06
#define STATUS_FP     (int) 0x07
#define STATUS_TN     (int) 0x08
#define STATUS_UNDEF  (int) 0x09

//const char * STATUS_STR[] = { "correct dest.", "wrong dest. > orig. server", "fp detect. > orig. server", "dropped"};

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "rid-router.h"

class Path_State {

    public:

        Path_State() {}
        Path_State(RID_Router * router, int request_size);
        ~Path_State() {}

        void set_in_fptree_prob(std::vector<__float080> * prob, int prob_size);
        void set_in_fptree_prob(int f, __float080 prob);
        void set_ingress_iface_prob(__float080 prob);

        std::vector<__float080> * get_in_fptree_prob();
        __float080 get_in_fptree_prob(uint8_t f);
        __float080 get_ingress_iface_prob();

        void set_eop();
        bool is_eop();

        void set_path_length(int length);
        void set_path_status(int status);
        void set_path_prob(__float080 prob);
        int get_path_length();
        int get_path_status();
        __float080 get_path_prob();

        void set_event(int event, __float080 prob);
        int get_event();
        __float080 get_event_prob();

        RID_Router * get_router();
        char * to_string();

        void set_tree_bitmask(std::vector<uint8_t> * tree_bitmask, int tree_bitmask_size);
        std::vector<uint8_t> * get_tree_bitmask();
        int get_tree_bitmask_size();

        void set_ttl(int ttl);
        int get_ttl();

    private:

        // a pointer to the RID router associated with this path state
        RID_Router * router;
        // is this the end of a path (EOP)?
        bool eop;
        // path status, path length (in hops) and path prob
        int path_length;
        int path_status;
        __float080 path_prob;
        // event & event prob
        int event;
        __float080 event_prob;
        // ingress probabilities : prob of a router receiving a request 
        // over a 'prefix tree' of size p
        std::vector<__float080> in_fptree_prob;
        __float080 ingress_iface_prob;
        // tree bitmask of the path
        std::vector<uint8_t> * tree_bitmask;
        int tree_bitmask_size;
        // ttl of packet
        int ttl;
};

#endif