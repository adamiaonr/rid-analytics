#ifndef PATH_STATE_HH
#define PATH_STATE_HH

#define MAX_PATH_STATE_STRING_SIZE 128

// possible outcomes
#define OUTCOME_UNDEF       0x00
#define OUTCOME_MULTI_HITS  0x01
#define OUTCOME_NO_HITS     0x02
#define OUTCOME_TP          0x04
#define OUTCOME_FP          0x08

//const char * OUTCOME_STR[] = { "correct dest.", "wrong dest. > orig. server", "fp detect. > orig. server", "dropped"};

#include <math.h>
#include <stdlib.h>
#include <stdio.h>

#include "rid-router.h"

class Path_State {

    public:

        Path_State() {}
        Path_State(__float080 ingress_prob, int path_length);
        ~Path_State() {}

        void set_ingress_prob(__float080 ingress_prob);
        __float080 get_ingress_prob();

        void set_eop();
        bool get_eop();

        void set_path_length(int path_length);
        int get_path_length();

        void set_outcome(uint8_t outcome);
        uint8_t get_outcome();

    private:

        // the probability of having a request reach this state
        __float080 ingress_prob;
        // is this the end of an RID packet path?
        bool is_eop;
        // length of the path, in hops
        int path_length;
        // outcome
        uint8_t outcome;
};

#endif