#include "path-state.h"

Path_State::Path_State(__float080 ingress_prob, int path_length) {

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
