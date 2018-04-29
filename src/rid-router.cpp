#include "rid-router.h"

const char * EVENT_STR[] = { 
    "EVENT_NLM", 
    "EVENT_MLM", 
    "EVENT_LLM", 
    "EVENT_SLM"};

const char * CUMULATIVE_PROB_MODE_STR[] = { 
    "MODE_EI_EXCLUSIVE", 
    "MODE_MIS", 
    "MODE_LI",
    "MODE_EI_INCLUSIVE"};

int RID_Router::init(
    std::string router_id,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t req_size,
    uint16_t bf_size,
    int mode) {

    if (this->initialized)
        return 0;

    this->id = router_id;

    this->req_size = req_size;
    this->bf_size = bf_size;
    // fwd table parameters
    std::cout << "RID_Router::init() : [INFO] iface_num = " << (int) iface_num << std::endl;
    this->iface_num = iface_num;
    this->fwd_table_size = fwd_table_size;
    // initialize forwarding table w/ 1 row per iface
    this->fwd_table = std::vector<RID_Router::fwd_table_row>(iface_num);
    // blocked ifaces initialized to false
    this->blocked_ifaces = std::vector<bool>(iface_num, false);
    // initialize iface on fp tree of size t array
    this->iface_on_fptree = std::vector< std::vector<bool> >(this->iface_num, std::vector<bool>(req_size, false));

    // initialize router's probability module
    this->prob_mod = new Prob((__float080) this->bf_size, (__float080) this->req_size, this->iface_num);

    // mark the router as initialized
    this->initialized = true;

    return 0;
}

RID_Router::RID_Router(
    std::string router_id,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t req_size,
    uint16_t bf_size,
    int mode) {

    this->initialized = false;
    this->init(router_id, fwd_table_size, iface_num, req_size, bf_size, mode);
}

RID_Router::~RID_Router() {}

int RID_Router::add_fwd_table_entry(
    int i, 
    __float080 iface_proportion, 
    std::map<int, __float080> * f_distr,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,
    RID_Router * next_hop_router,
    int next_hop_iface) {

    if ((i < 0) || (i > (this->iface_num - 1))) {
        std::cerr << "RID_Router::add_fwd_table_entry() : [ERROR] invalid iface index: "
            << i << " vs. " << (int) this->iface_num << std::endl;
        return -1;
    }

    this->fwd_table[i].num_entries = (__float080) this->fwd_table_size * iface_proportion;
    std::cout << "RID_Router::add_fwd_table_entry() : [INFO] iface[" << i 
        << "].num_entries = " << this->fwd_table[i].num_entries << std::endl;
    
    this->fwd_table[i].iface = i;
    // iface_proportion : % of table entries associated w/ this iface
    this->fwd_table[i].iface_proportion = iface_proportion;
    // distribution of sizes among the entries associated w/ this iface
    this->fwd_table[i].f_distr = std::vector<__float080>(this->req_size, 0.0);

    for (int f = 0; f < this->req_size; f++) {
        this->fwd_table[i].f_distr[f] = (*f_distr)[f + 1];
        // set bit in the bitmask of fwd entry sizes if f_distr[f] > 0.0
        if (this->fwd_table[i].f_distr[f] > 0.0)
            this->fwd_table[i].f_bitmask |= (1 << f);
    }

    // copy the map of tree bitmasks into the router's iface
    this->fwd_table[i].tree_bitmasks = (*tree_bitmasks);

    // next hop information: remote router pointer and iface
    this->fwd_table[i].next_hop.router = next_hop_router;
    this->fwd_table[i].next_hop.iface = next_hop_iface;

    return 0;
}

Prob::fp_data * RID_Router::get_fp_data(
    RID_Router * router, 
    uint8_t i,
    bool iface_complement) {

    Prob::fp_data * iface_fp_data = new Prob::fp_data();

    if (!iface_complement) {

        iface_fp_data->tp_size = router->tp_sizes[i];
        iface_fp_data->num_entries = (__float080) router->fwd_table[i].num_entries;
        iface_fp_data->f_distr = router->fwd_table[i].f_distr;
        iface_fp_data->f_r_distr = router->fwd_table[i].f_r_distr;
        iface_fp_data->entry_prop = router->fwd_table[i].iface_proportion;
        iface_fp_data->tree_bitmasks = router->fwd_table[i].tree_bitmasks;
        iface_fp_data->on_fptree = router->iface_on_fptree[i];
        iface_fp_data->is_blocked = router->blocked_ifaces[i];

    } else {

        iface_fp_data->f_distr = std::vector<__float080>(this->req_size, 0.0);
        iface_fp_data->f_r_distr = std::vector<__float080>(this->req_size, 0.0);
        iface_fp_data->tree_bitmasks = router->fwd_table[i].tree_bitmasks;
        // clear iface fp data's tree bitmask
        std::map<int, std::vector<uint8_t> >::iterator itr;
        for (itr = iface_fp_data->tree_bitmasks.begin(); itr != iface_fp_data->tree_bitmasks.end(); itr++)
            std::fill(itr->second.begin(), itr->second.end(), 0x00);

        // sum all the fp data, for all ifaces other than i
        __float080 num_valid_ifaces = 0.0;
        for (uint8_t k = 1; k < router->iface_num; k++) {

            // don't consider iface i, or blocked ifaces
            if ((k == i) || (router->fwd_table[k].num_entries < 1)) continue;
            if (router->blocked_ifaces[k]) continue;

            // the tp is the max of the ifaces other than i
            iface_fp_data->tp_size = std::max(iface_fp_data->tp_size, router->tp_sizes[k]);
            // keep adding up the nr of entries
            iface_fp_data->num_entries += (__float080) router->fwd_table[k].num_entries;

            // f_distr and f_r_distr vectors
            // use std::transform() to get the elementwise sum of 2 vectors 
            // (result saved in first vector)
            std::transform(
                iface_fp_data->f_distr.begin(), iface_fp_data->f_distr.end(), 
                router->fwd_table[k].f_distr.begin(), 
                iface_fp_data->f_distr.begin(), 
                std::plus<__float080>());

            std::transform(
                iface_fp_data->f_r_distr.begin(), iface_fp_data->f_r_distr.end(), 
                router->fwd_table[k].f_r_distr.begin(), 
                iface_fp_data->f_r_distr.begin(), 
                std::plus<__float080>());

            // fill the bitmask of ~i in the old fashioned way...
            for (itr = iface_fp_data->tree_bitmasks.begin(); itr != iface_fp_data->tree_bitmasks.end(); itr++) {
                int tb_size = itr->first;
                for (unsigned int j = 0; j < router->fwd_table[k].tree_bitmasks[tb_size].size(); j++)
                    itr->second[j] |= router->fwd_table[k].tree_bitmasks[tb_size][j];
            }

            iface_fp_data->entry_prop += router->fwd_table[k].iface_proportion;

            iface_fp_data->on_fptree = std::vector<bool>(this->req_size, false);
            for (unsigned int j = 0; j < this->req_size; j++) {
                // std::cout << "RID_Router::get_fp_data() : [INFO] " << 
                //     "\n\ton_fptree[~" << (int) i << "][" << j << "] : " << iface_fp_data->on_fptree[j] << 
                //     "\n\ton_fptree[" << (int) k << "][" << j << "] : " << router->iface_on_fptree[k][j] << std::endl;
                iface_fp_data->on_fptree[j] = (iface_fp_data->on_fptree[j] || router->iface_on_fptree[k][j]);
                // std::cout << "RID_Router::get_fp_data() : [INFO] " << 
                //     "\n\ton_fptree[~" << (int) i << "][" << j << "] : " << iface_fp_data->on_fptree[j] << std::endl;
            }

            num_valid_ifaces++;
        }

        // normalize the iface_fp_data->f_distr and iface_fp_data->f_r_distr vectors
        for (uint8_t k; k < iface_fp_data->f_distr.size(); k++)
            iface_fp_data->f_distr[k] /= num_valid_ifaces;
        for (uint8_t k; k < iface_fp_data->f_r_distr.size(); k++)
            iface_fp_data->f_r_distr[k] /= num_valid_ifaces;
    }

    return iface_fp_data;
}

int RID_Router::forward(
    uint8_t ingress_iface,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,
    __float080 in_prob,
    std::vector<__float080> * in_fptree_prob,
    std::vector<std::vector<__float080> > & iface_probs,
    std::vector<std::vector<__float080> > & out_fptree_probs) {

    std::cout << "RID_Router::forward() : [INFO] *** forwarding *** " << 
        "\n\tfrom : " << this->get_next_hop(ingress_iface).router->get_id() << 
        "\n\tinto : " << this->get_id() << "[" << (int) ingress_iface << "]" << std::endl;

    // check which ifaces are eligible to be on prefix trees (of any size)
    for (uint8_t i = 0; i < this->iface_num; i++) {
        std::fill(this->iface_on_fptree[i].begin(), this->iface_on_fptree[i].end(), false);
        is_iface_on_fptree(i, tree_bitmasks);
    }

    // mark blocked ifaces for basic loop prevention
    // don't forward over ifaces w/ no entries associated w/ them
    for (uint8_t i = 0; i < this->iface_num; i++) {

        this->blocked_ifaces[i] = false;
        
        if (this->fwd_table[i].num_entries == 0)
            this->blocked_ifaces[i] = true;
    }
    // don't foward over ingress iface
    this->blocked_ifaces[ingress_iface] = true;

    std::cout << "RID_Router::forward() : [INFO] blocked ifaces {";
    for (uint8_t k = 0; k < this->blocked_ifaces.size(); k++) {
        std::cout << "[" << (int) k << "]:" << (int) this->blocked_ifaces[k] << ", ";        
    }
    std::cout << "}" << std::endl;

    // save the ingress probabilities for each diff. prefix tree size (i.e. 
    // the probability that a packet is already bound to a tree of size p)

    std::cout << "RID_Router::forward() : P(in) = " << in_prob << std::endl;
    std::cout << "RID_Router::forward() : P(in_fptree = p) :" << std::endl; 
    for (uint8_t fptree_size = 0; fptree_size <= this->req_size; fptree_size++)
        std::cout << "\tP(p = " << (int) fptree_size << "]) = " << (*in_fptree_prob)[fptree_size] << std::endl;

    // get the iface probs P(I = i), i.e. the probability of having iface i 
    // chosen as an egress iface 

    // as input, we pass: 
    //  - iface data from the fwd_table[] which has an influence in fp rate 
    //    and/or iface prob calculation
    //  - in_prob : probability of having a request forwarded 
    //              into this router by the previous router.
    //  - in_fptree_prob : probabilities of having the request coming into this 
    //                     router while 'bound' to a fp tree of size p
    std::vector<std::vector<Prob::fp_data *> > iface_fp_data(2);
    for (uint8_t i = 0; i < this->iface_num; i++) {

        iface_fp_data[0].push_back(get_fp_data(this, i));
        iface_fp_data[1].push_back(get_fp_data(this, i, true));
    }

    if (this->prob_mod->calc_probs(
        &(iface_fp_data),
        tree_bitmasks,
        in_prob,
        in_fptree_prob, 
        iface_probs, 
        out_fptree_probs) < 0) return -1;

    std::cout << "RID_Router::forward() : [INFO] iface probs :" << std::endl;
    for(uint8_t i = 0; i < this->iface_num; i++)
        std::cout << "\tP(I = " << (int) i << ") = " << iface_probs[i][0] << " : " << iface_probs[i][1] << " : " << iface_probs[i][2] << std::endl;

    return 0;
}

// objective: find if iface i is part of a fp tree which started at a 
//            previous router.
//
// how do we encode the fp trees? 
// basically, using a bitmask. 
// the k-th bit in the bitmask tells if the node w/ id equal to k can be 
// reached by the current router or iface. 
// e.g. if the bitmask is '10001000', then only the 
// nodes w/ ids 3 and 7 can be reached from the current router.
//
// there are 2 bitmasks in this context:
//  - tree bitmask      : encodes nodes reachable by current router 
//  - iface[i] bitmask  : encodes nodes reachable by iface i
//
// algorithm:
//  - for every bit k in iface[i] bitmask
//      - if iface[i][k] == tree_bitmask[k], return true
//          - i.e. iface i can be used to reach a node in the tree bitmask
//  - return false : iface i can't reach any of the nodes in the tree bitmask
//
// why is this useful? 
//
void RID_Router::is_iface_on_fptree(
    int i, 
    std::map<int, std::vector<uint8_t> > * tree_bitmasks) {

    // std::cout << "RID_Router::is_iface_on_fptree() : [INFO] is " << i << " on prefix tree?";

    std::map<int, std::vector<uint8_t> >::iterator itr;
    for (itr = (*tree_bitmasks).begin(); itr != (*tree_bitmasks).end(); itr++) {

        int t = itr->first;
        if (t == 0) continue;
    
        for (uint8_t k = 0; k < this->fwd_table[i].tree_bitmasks[t].size(); k++) {

            uint8_t iface_tree_bitmask = this->fwd_table[i].tree_bitmasks[t][k];
            // std::cout << "\n\t tree_bitmask vs. iface_tree_bitmask : " 
            //     << (int) (*tree_bitmasks)[t][k] << " <-> " << (int) iface_tree_bitmask
            //     << std::endl;

            if ((*tree_bitmasks)[t][k] & iface_tree_bitmask) {
                // std::cout << "\n\t IT IS!" << std::endl;
                this->iface_on_fptree[i][t] = true;
            }
        }
    }
}