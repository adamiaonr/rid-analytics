#ifndef PATH_STATE_HH
#define PATH_STATE_HH

#define MAX_PATH_STATE_STRING_SIZE 128

// possible outcomes
#define OUTCOME_CORRECT_DELIVERY    0x00
#define OUTCOME_INCORRECT_DELIVERY  0x01
#define OUTCOME_FALLBACK_DELIVERY   0x02
#define OUTCOME_FALLBACK_RELAY      0x03
#define OUTCOME_INTERMEDIATE_TP     0x04
#define OUTCOME_INTERMEDIATE_FP     0x05
#define OUTCOME_UNDEF               0x06

//const char * OUTCOME_STR[] = { "correct dest.", "wrong dest. > orig. server", "fp detect. > orig. server", "dropped"};

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "rid-router.h"

class Path_State {

    public:

        Path_State() {}
        Path_State(RID_Router * rid_router, int request_size);
        ~Path_State() {}

        void set_final_prob(__float080 prob);
        __float080 get_final_prob();

        void set_ingress_ptree_prob(__float080 * prob, int prob_size);
        void set_ingress_ptree_prob(int f, __float080 prob);
        void set_ingress_iface_prob(__float080 prob);
        
        __float080 * get_ingress_ptree_prob();
        __float080 get_ingress_ptree_prob(uint8_t f);
        __float080 get_ingress_iface_prob();

        void set_eop();
        bool get_eop();

        void set_path_length(int path_length);
        int get_path_length();

        void set_outcome(uint8_t outcome);
        uint8_t get_outcome();

        RID_Router * get_router();

        char * to_string();

    private:

        // a pointer to the RID router associated with this path state
        RID_Router * rid_router;

        // request sizes used for this run
        int request_size;

        // is this the end of an RID packet path?
        bool is_eop;

        // length of the path, in hops
        int path_length;

        // outcome
        uint8_t outcome;

        // single probability for the state, i.e. interface event. in the 
        // case of a EI event, this probability is the sum of the ingress 
        // 'prefix tree' probabilities for all sizes
        __float080 final_prob;

        // ingress probabilities : prob of a router receiving a request 
        // over a 'prefix tree' of size p
        __float080 * ingress_ptree_prob;
        __float080 ingress_iface_prob;
};

#endif