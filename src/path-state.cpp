#include "path-state.h"

Path_State::Path_State(
    RID_Router * rid_router, 
    __float080 ingress_prob, 
    int path_length) {

    this->rid_router = rid_router;
    this->ingress_prob = ingress_prob;
    this->path_length = path_length;
    this->is_eop = false;
    this->outcome = OUTCOME_UNDEF;
}

void Path_State::set_ingress_prob(__float080 ingress_prob) {
    this->ingress_prob = ingress_prob;
}

__float080 Path_State::get_ingress_prob() {
    return this->ingress_prob;
}

void Path_State::set_eop() { 
    this->is_eop = true; 
}

bool Path_State::get_eop() { 
    return this->is_eop; 
}

void Path_State::set_path_length(int path_length) { 
    this->path_length = path_length; 
}

int Path_State::get_path_length() { 
    return this->path_length; 
}

void Path_State::set_outcome(uint8_t outcome) { 
    this->outcome = outcome; 
}

uint8_t Path_State::get_outcome() { 
    return this->outcome; 
}

RID_Router * Path_State::get_router() { 
    return this->rid_router; 
}

char * Path_State::to_string() {

    // using calloc() since node_str will be returned
    char * node_str = (char *) calloc(MAX_PATH_STATE_STRING_SIZE, sizeof(char));
    char _outcome[MAX_PATH_STATE_STRING_SIZE];

    switch(this->outcome) {

        case OUTCOME_MULTI_HITS:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "mult. iface hits");
            break;
        case OUTCOME_NO_HITS:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "no hits");
            break;
        case OUTCOME_TP:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "tp hit");
            break;
        case OUTCOME_FP:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "fp hit");
            break;
        default:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "unknown");
            break;
    }

    snprintf(
        node_str, MAX_PATH_STATE_STRING_SIZE, 
        "ROUTER[%d][%d] : [PROB : %-.5LE][PATH_LENGTH : %-5d][OUTCOME: %s]",
        this->rid_router->get_height(), this->rid_router->get_width(),
        this->ingress_prob, this->path_length, _outcome);

    return node_str;
}
