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
    int mm_mode) {

    if (this->initialized)
        return 0;

    this->id = router_id;

    this->req_size = req_size;
    this->bf_size = bf_size;

    // fwd table parameters
    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    // initialize forwarding table w/ 1 row per iface
    this->fwd_table = std::vector<RID_Router::fwd_table_row>(iface_num);

    // // iface event probabilities. events are: 
    // //  - EVENT_NLM: no link matches
    // //  - EVENT_MLM: multiple link matches
    // //  - EVENT_LLM: local link match
    // //  - EVENT_SLM: single link match (other than local)
    // this->iface_events_prob = (__float080 *) calloc(EVENT_NUM, sizeof(__float080));

    // // initialize the in_fptree_prob[p] array. holds the prob of having 
    // // request entering a router on a prefix tree of size p.
    // this->in_fptree_prob = (__float080 *) calloc(this->req_size + 1, sizeof(__float080));
    // // initialize the egress_ptree_prob[i][p] array. holds the prob of having 
    // // request leaving a router on iface i, on a prefix tree of size p.
    // this->egress_ptree_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    // for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
    //     this->egress_ptree_prob[_iface] = (__float080 *) calloc(this->req_size + 1, sizeof(__float080));

    // // initialize the egress_iface_prob[i][f] array. holds the probability 
    // // of having a request leaving a router over iface i, due to a match of 
    // // size f.
    // this->egress_iface_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    // for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
    //     this->egress_iface_prob[_iface] = (__float080 *) calloc(this->req_size + 1, sizeof(__float080));

    this->mm_mode = mm_mode;

    // initialize router's probability module
    this->prob_mod = new Prob((__float080) this->bf_size, (__float080) this->req_size, this->iface_num);

    this->initialized = true;

    return 0;
}

RID_Router::RID_Router(
    std::string router_id,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t req_size,
    uint16_t bf_size,
    int mm_mode) {

    this->init(router_id, fwd_table_size, iface_num, req_size, bf_size, mm_mode);
}

RID_Router::~RID_Router() {}

int RID_Router::add_fwd_table_entry(
    int i, 
    __float080 iface_proportion, 
    std::map<int, __float080> * f_distr,
    std::vector<uint8_t> * tree_bitmask,
    RID_Router * next_hop_router,
    int next_hop_iface) {

    if ((i < 0) || (i > (this->iface_num - 1))) {
        std::cerr << "RID_Router::add_fwd_table_entry() : [ERROR] invalid iface index: "
            << i << std::endl;
        return -1;
    }

    this->fwd_table[i].num_entries = (__float080) this->fwd_table_size * iface_proportion;
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

    // allocate memory for and fill the tree bitmask, containing the 
    // source trees included in this iface
    this->fwd_table[i].tree_bitmask = *tree_bitmask;

    // next hop information: remote router pointer and iface
    this->fwd_table[i].next_hop.router = next_hop_router;
    this->fwd_table[i].next_hop.iface = next_hop_iface;

    return 0;
}

Prob::fp_data RID_Router::get_fp_data(
    RID_Router * router, 
    uint8_t i, 
    bool anti) {

    Prob::fp_data iface_fp_data;

    if (!anti) {

        iface_fp_data.tp_size = router->tp_sizes[i];
        iface_fp_data.num_entries = (__float080) router->fwd_table[i].num_entries;
        // pointers to f distr.
        iface_fp_data.f_distr = router->fwd_table[i].f_distr;
        iface_fp_data.f_r_distr = router->fwd_table[i].f_r_distr;

        iface_fp_data.on_fptree = router->iface_on_fptree[i];

    } else {

        // sum all the fp data, for all ifaces other than i
        for (uint8_t k = 0; k < router->iface_num; k++) {

            if ((k == i) || (router->fwd_table[k].num_entries < 1)) continue;

            // the tp is the max of the ifaces other than i
            iface_fp_data.tp_size = std::max(iface_fp_data.tp_size, router->tp_sizes[k]);
            // keep adding up the nr of entries
            iface_fp_data.num_entries += (__float080) router->fwd_table[k].num_entries;
            // elementwise sum of 2 vectors (result saved in first vector)
            std::transform(
                iface_fp_data.f_distr.begin(), iface_fp_data.f_distr.end(), 
                router->fwd_table[k].f_r_distr.begin(), 
                iface_fp_data.f_distr.begin(), 
                std::plus<__float080>());

            std::transform(
                iface_fp_data.f_r_distr.begin(), iface_fp_data.f_r_distr.end(), 
                router->fwd_table[k].f_r_distr.begin(), 
                iface_fp_data.f_r_distr.begin(), 
                std::plus<__float080>());

            iface_fp_data.on_fptree |= router->iface_on_fptree[k];
        }
    }

    return iface_fp_data;
}

int RID_Router::forward(
    uint8_t ingress_iface,
    std::vector<uint8_t> * tree_bitmask,
    __float080 ingress_prob,
    __float080 * in_fptree_prob) {

    std::cout << "RID_Router::forward() : [INFO] forwarding " << 
        "\n\tfrom : " << this->get_next_hop(ingress_iface).router->get_id() << 
        "\n\tinto : " << this->get_id() << "[" << (int) ingress_iface << "]" << std::endl;

    // check which ifaces are eligible to be on prefix trees (of any size)
    this->iface_on_fptree.clear();
    for (uint8_t i = 0; i < this->iface_num; i++)
        this->iface_on_fptree.push_back(is_iface_on_fptree(i, tree_bitmask));

    // // pre-calculate fptree iface probabilities
    // calc_iface_on_fptree_probs();

    // basic loop prevention : keep a list of ifaces over which the request
    // should not be forwarded
    this->blocked_ifaces.clear();
    // don't foward over ingress iface
    this->blocked_ifaces.insert(ingress_iface);
    // don't forward over ifaces w/ no entries associated w/ them
    for (uint8_t i = 0; i < this->iface_num; i++)
        if (this->fwd_table[i].num_entries == 0)
            this->blocked_ifaces.insert(i);

    std::cout << "RID_Router::forward() : [INFO] not forwarding over ifaces {";
    for (std::set<uint8_t>::iterator it = this->blocked_ifaces.begin();
        it != this->blocked_ifaces.end();
        ++it) {

        std::cout << (int) (*it) << ", ";        
    }
    std::cout << "}" << std::endl;

    // save the ingress probabilities for each diff. prefix tree size (i.e. 
    // the probability that a packet is already bound to a tree of size p)

    // ingress probabilities, what are these?
    //  - ingress_prob : probability of having the request forwarded 
    //                   INTO this router by the previous router, regardless 
    //                   of the cause.
    //
    //  - in_fptree_prob : probability of having the request coming in while 
    //                         'stuck' in a 'false positive' prefix tree of 
    //                         size p, with 1 <= p <= req_size.
    // this->ingress_prob = ingress_prob;
    // std::cout << "RID_Router::forward() : P(INGRESS) = " << this->ingress_prob << std::endl;

    // for (uint8_t fptree_size = 0; fptree_size <= this->req_size; fptree_size++) {
    //     this->in_fptree_prob[fptree_size] = in_fptree_prob[fptree_size];
    //     std::cout << "RID_Router::forward() : P(INGRESS_PTREE[" 
    //         << (int) ptree_size << "]) = " << this->in_fptree_prob[ptree_size] << std::endl;
    // }
    
    // 1) calculate largest match probabilities
    std::vector<std::vector<Prob::fp_data> > iface_fp_data(2, std::vector<Prob::fp_data> (this->iface_num));
    for (uint8_t i = 0; i < this->iface_num; i++) {
        iface_fp_data[0].push_back(get_fp_data(this, i));
        iface_fp_data[1].push_back(get_fp_data(this, i, true));
    }

    if (this->prob_mod->calc_lm_prob(&iface_fp_data) < 0) return -1;

    // 1.2) print egress_ptree_prob[i][p], i.e. the 
    // probability of having the request bound to a fp prefix tree of size p, 
    // over iface i

    // 2) calculate the joint probability distribution of the interface largest 
    // matches, i.e. the probability of a forwarding event. a forwarding event 
    // is represented by a 1 x iface_num vector, e.g. [1, 2, 5], which means 
    // 'largest match of size 1 at iface 0, lm of size 2 at iface 1 and 
    // lm of 5 at iface 3'. this allows us to answer several questions about 
    // forwarding events:
    //  - how likely is it to forward the request over iface i?
    //  - how likely is it to forward the request over multiple ifaces?
    //  - how likely is it to have no matches at all, i.e. event {0, 0, 0, 0}?

    // 2.1) initializations

    // 2.1.1) initialize matrix on which to save the joint largest match 
    // probabilities
    // FIXME: this is incredibly complex...

    // 2.1.2) initialize iface_events_prob

    // 2.1.3) initialize egress iface probs

    // 2.2) compute the fwd event probs given that the request is bound to 
    // a fp tree of size s. under such conditions, the possible events are 
    // restricted to:
    //  - fp tree of size s continues from the previous router, OR (nor XOR)
    //  - new fp tree, of size >= s, starts at this router

    return 0;
}

bool RID_Router::is_iface_on_fptree(int i, std::vector<uint8_t> * tree_bitmask) {

    std::cout << "RID_Router::is_iface_on_fptree() : [INFO] is " << i 
        << " on prefix tree?";

    for (uint8_t t = 0; t < this->fwd_table[i].tree_bitmask.size(); t++) {

        uint8_t iface_tree_bitmask = this->fwd_table[i].tree_bitmask[t];

        std::cout << "\n\t tree_bitmask vs. iface_tree_bitmask : " 
            << (int) (*tree_bitmask)[t] << " <-> " << (int) iface_tree_bitmask;

        if ((*tree_bitmask)[t] & iface_tree_bitmask) {
            std::cout << "\n\t IT IS!" << std::endl;
            return true;
        }
    }

    std::cout << "\n\t IT'S NOT..." << std::endl;
    return false;
}