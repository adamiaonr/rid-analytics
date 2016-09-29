#include "path-state.h"

Path_State::Path_State(RID_Router * rid_router, int request_size) {

    this->rid_router = rid_router;
    this->request_size = request_size;
    // requests can come in by hitting FP matches of diff. sizes at previous 
    // routers
    this->final_prob = 0.0;
    this->ingress_ptree_prob = (__float080 *) calloc(request_size + 1, sizeof(__float080));
    this->is_eop = false;
    this->outcome = OUTCOME_UNDEF;
}

void Path_State::set_final_prob(__float080 prob) {

    this->final_prob = prob;
}

__float080 Path_State::get_final_prob() {

    return this->final_prob;
}

void Path_State::set_ingress_ptree_prob(__float080 * prob, int prob_size) {

    for (int f = 0; f < prob_size + 1; f++)
        this->ingress_ptree_prob[f] = prob[f];
}

void Path_State::set_ingress_ptree_prob(int f, __float080 prob) {

    this->ingress_ptree_prob[f] = prob;
}

void Path_State::set_ingress_iface_prob(__float080 prob) {

    this->ingress_iface_prob = prob;
}

__float080 * Path_State::get_ingress_ptree_prob() {

    return this->ingress_ptree_prob;
}

__float080 Path_State::get_ingress_ptree_prob(uint8_t f) {

    return this->ingress_ptree_prob[f];
}

__float080 Path_State::get_ingress_iface_prob() {

    return this->ingress_iface_prob;
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
    char * _node_str = (char *) calloc(MAX_PATH_STATE_STRING_SIZE, sizeof(char));
    char _outcome[MAX_PATH_STATE_STRING_SIZE];

    switch(this->outcome) {

        case OUTCOME_CORRECT_DELIVERY:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "CORRECT_DELIVERY");
            break;
        case OUTCOME_INCORRECT_DELIVERY:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "INCORRECT_DELIVERY");
            break;
        case OUTCOME_FALLBACK_DELIVERY:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "FALLBACK_DELIVERY");
            break;
        case OUTCOME_FALLBACK_RELAY:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "FALLBACK_RELAY");
            break;
        case OUTCOME_INTERMEDIATE_TP:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "INTERMEDIATE_TP");
            break;
        case OUTCOME_INTERMEDIATE_FP:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "INTERMEDIATE_FP");
            break;
        default:
            snprintf(_outcome, MAX_PATH_STATE_STRING_SIZE, "UNKNOWN");
            break;
    }

    snprintf(
        _node_str, MAX_PATH_STATE_STRING_SIZE, 
        "ROUTER[%d][%d] : [PROB : %-.5LE][PATH_LENGTH : %-5d][OUTCOME: %s]",
        this->rid_router->get_height(), this->rid_router->get_width(),
        this->final_prob, this->path_length, _outcome);

    return _node_str;
}
