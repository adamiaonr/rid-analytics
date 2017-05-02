#include "path-state.h"

Path_State::Path_State(RID_Router * router, int request_size) {

    this->router = router;
    // path info
    this->path_length = 0;
    this->path_status = (int) OUTCOME_UNDEF;
    this->path_prob = 0.0;
    // event & event prob
    this->event = (int) EVENT_UNKNOWN;
    this->event_prob = 0.0;

    this->ingress_ptree_prob = (__float080 *) calloc(request_size + 1, sizeof(__float080));
    this->eop = false;
}

void Path_State::set_path_length(int length) {
    this->path_length = length;
}

int Path_State::get_path_length() {
    return this->path_length;
}

void Path_State::set_path_status(int status) {
    this->path_status = status;
}

int Path_State::get_path_status() {
    return this->path_status;
}

void Path_State::set_path_prob(__float080 prob) {
    this->path_prob = prob;
}

__float080 Path_State::get_path_prob() {
    return this->path_prob;
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
    this->eop = true; 
}

bool Path_State::is_eop() { 
    return this->eop; 
}

void Path_State::set_event(int event, __float080 prob) { 
    this->event = event; 
    this->event_prob = prob;
}

int Path_State::get_event() { 
    return this->event; 
}

__float080 Path_State::get_event_prob() { 
    return this->event_prob; 
}

void Path_State::set_tree_bitmask(uint8_t * tree_bitmask, int tree_bitmask_size) {
    this->tree_bitmask_size = tree_bitmask_size;
    this->tree_bitmask = tree_bitmask;
}

int Path_State::get_tree_bitmask_size() { 
    return this->tree_bitmask_size; 
}

uint8_t * Path_State::get_tree_bitmask() { 
    return this->tree_bitmask; 
}

void Path_State::set_ttl(int ttl) { 
    this->ttl = ttl;
}

int Path_State::get_ttl() { 
    return this->ttl; 
}

RID_Router * Path_State::get_router() { 
    return this->router; 
}

char * Path_State::to_string() {

    // using calloc() since node_str will be returned
    char * node_str = (char *) calloc(MAX_PATH_STATE_STRING_SIZE, sizeof(char));
    char status[MAX_PATH_STATE_STRING_SIZE];

    switch(this->path_status) {

        case OUTCOME_CORRECT_DELIVERY:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "CORRECT_DELIVERY");
            break;
        case OUTCOME_INCORRECT_DELIVERY:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "INCORRECT_DELIVERY");
            break;
        case OUTCOME_FALLBACK_DELIVERY:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "FALLBACK_DELIVERY");
            break;
        case OUTCOME_FALLBACK_RELAY:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "FALLBACK_RELAY");
            break;
        case OUTCOME_PACKET_DROP:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "PACKET_DROP");
            break;
        case OUTCOME_TTL_DROP:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "TTL_DROP");
            break;
        case STATUS_TP:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "STATUS_TP");
            break;
        case STATUS_FP:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "STATUS_FP");
            break;
        case STATUS_TN:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "STATUS_TN");
            break;
        default:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "UNKNOWN");
            break;
    }

    snprintf(
        node_str, MAX_PATH_STATE_STRING_SIZE, 
        "ROUTER[%s] : [PROB : %-.5LE][PATH_LENGTH : %-5d][OUTCOME: %s]",
        this->router->get_id().c_str(),
        this->path_prob, this->path_length, status);

    return node_str;
}
