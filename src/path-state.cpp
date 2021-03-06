#include "path-state.h"

Path_State::Path_State(RID_Router * router, int req_size) {

    this->router = router;
    // path info
    this->path_length = 0;
    this->path_status = (int) OUTCOME_UNDEF;
    this->path_prob = 0.0;
    // event & event prob
    this->event = (int) 0x07;
    this->event_prob = 0.0;

    this->in_fptree_prob = std::vector<__float080>(req_size + 1, 0.0);

    this->eop = false;
}

void Path_State::set_length(int length) {
    this->path_length = length;
}

int Path_State::get_length() {
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

void Path_State::set_in_fptree_prob(std::vector<__float080> * prob, uint8_t prob_size) {
    for (unsigned int f = 0; f < (*prob).size(); f++)
        this->in_fptree_prob[f] = (*prob)[f];
}

void Path_State::set_in_fptree_prob(int f, __float080 prob) {
    this->in_fptree_prob[f] = prob;
}

void Path_State::set_in_iface_prob(__float080 prob) {
    this->in_iface_prob = prob;
}

std::vector<__float080> * Path_State::get_in_fptree_prob() {
    return &(this->in_fptree_prob);
}

__float080 Path_State::get_in_fptree_prob(uint8_t f) {
    return this->in_fptree_prob[f];
}

__float080 Path_State::get_in_iface_prob() {

    return this->in_iface_prob;
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

void Path_State::set_tree_bitmasks(
    std::map<int, std::vector<uint8_t> > * tree_bitmasks) {

    this->tree_bitmasks = tree_bitmasks;
}

std::map<int, std::vector<uint8_t> > * Path_State::get_tree_bitmasks() { 
    return this->tree_bitmasks; 
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
        case OUTCOME_HARD_FALLBACK_DELIVERY:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "HARD_FALLBACK_DELIVERY");
        case OUTCOME_SOFT_FALLBACK_FWD:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "SOFT_FALLBACK_FWD");
            break;
        case OUTCOME_FEEDBACK_DELIVERY:
            snprintf(status, MAX_PATH_STATE_STRING_SIZE, "FEEDBACK_DELIVERY");
            break;
        // case OUTCOME_FALLBACK_RELAY:
        //     snprintf(status, MAX_PATH_STATE_STRING_SIZE, "FALLBACK_RELAY");
        //     break;
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
