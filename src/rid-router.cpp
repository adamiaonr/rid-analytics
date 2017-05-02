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

RID_Router::RID_Router(
    std::string router_id,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max,
    uint16_t bf_size,
    int mm_mode) {

    // router id
    this->id = router_id;
    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    this->f_max = f_max;
    this->bf_size = bf_size;

    // initialize forwarding table
    this->fwd_table = (RID_Router::fwd_table_row *) calloc(iface_num, sizeof(RID_Router::fwd_table_row)); 
    // to indicate an empty entry in the table, set iface to -1
    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
        this->fwd_table[i].iface = -1;

    // lpm_pmf is a 2D array
    this->lpm_pmf = (RID_Router::lpm_pmf_row **) calloc(this->iface_num, sizeof(RID_Router::lpm_pmf_row));

    // initialize the iface_events_pmf[] array
    this->iface_events_pmf = (__float080 *) calloc(EVENT_NUM, sizeof(__float080));

    // initialize the ingress_ptree_prob[] array
    this->ingress_ptree_prob = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    // initialize the egress_ptree_prob[iface][] array
    this->egress_ptree_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_ptree_prob[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    // initialize the egress_iface_prob[iface][] array
    this->egress_iface_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_iface_prob[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    this->initialized = true;
    this->mm_mode = mm_mode;
}

RID_Router::~RID_Router() {

    // free forwarding table
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++) {

        free(this->fwd_table[_iface].f_distribution);
    }

    free(this->fwd_table);

    // free lpm_pmf
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++) {

        if (this->lpm_pmf[_iface] == NULL)
            continue;

        for (uint8_t _f = 0; _f < this->f_max; _f++) {

            free(this->lpm_pmf[_iface][_f].lpm_pmf_prob);
        }

        free(this->lpm_pmf[_iface]);
    }

    free(this->iface_events_pmf);

    for (uint8_t _iface = IFACE_LOCAL + 1; _iface < this->iface_num; _iface++)
        free(this->egress_ptree_prob[_iface]);

    for (uint8_t _iface = IFACE_LOCAL + 1; _iface < this->iface_num; _iface++)
        free(this->egress_iface_prob[_iface]);
}

int RID_Router::init(
    std::string router_id,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max,
    uint16_t bf_size,
    int mm_mode) {

    if (this->initialized)
        return 0;

    // router id
    this->id = router_id;
    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    this->f_max = f_max;
    this->bf_size = bf_size;

    // initialize forwarding table
    this->fwd_table = (RID_Router::fwd_table_row *) calloc(iface_num, sizeof(RID_Router::fwd_table_row)); 
    // to indicate an empty entry in the table, set iface to -1
    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
        this->fwd_table[i].iface = -1;

    // lpm_pmf is a 2D array
    this->lpm_pmf = (RID_Router::lpm_pmf_row **) calloc(this->iface_num, sizeof(RID_Router::lpm_pmf_row));

    // initialize the iface_events_pmf[] array
    this->iface_events_pmf = (__float080 *) calloc(EVENT_NUM, sizeof(__float080));

    // initialize the ingress_ptree_prob[] array
    this->ingress_ptree_prob = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    // initialize the egress_ptree_prob[iface][] array
    this->egress_ptree_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_ptree_prob[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    // initialize the egress_iface_prob[iface][] array
    this->egress_iface_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_iface_prob[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    this->initialized = true;
    this->mm_mode = mm_mode;
    return 0;
}

int RID_Router::add_fwd_table_entry(
    int iface, 
    __float080 iface_proportion, 
    std::map<int, __float080> * size_dist,
    std::vector<uint8_t> * tree_bitmask,
    RID_Router * next_hop_router,
    int next_hop_iface) {

    if (iface < (int) IFACE_LOCAL || iface > (this->iface_num - 1)) {
        std::cerr << "RID_Router::add_fwd_table_entry() : [ERROR] invalid iface index: "
            << iface << std::endl;
        return -1;
    }

    this->fwd_table[iface].num_entries = (__float080) this->fwd_table_size * iface_proportion;
    this->fwd_table[iface].iface = iface;
    // iface_proportion : % of table entries associated w/ this iface
    this->fwd_table[iface].iface_proportion = iface_proportion;

    // distribution of sizes among the entries associated w/ this iface
    this->fwd_table[iface].f_distribution = 
        (__float080 *) calloc(this->f_max, sizeof(__float080));
    for (int f = 0; f < this->f_max; f++) {
        this->fwd_table[iface].f_distribution[f] = (*size_dist)[f + 1];
        // set bit in the bitmask of fwd entry sizes if f_distribution[f] > 0.0
        if (this->fwd_table[iface].f_distribution[f] > 0.0)
            this->fwd_table[iface].f_bitmask |= (1 << f);
    }

    // allocate memory for and fill the tree bitmask, containing the 
    // source trees included in this iface
    this->fwd_table[iface].tree_bitmask = 
        (uint8_t *) calloc((*tree_bitmask).size(), sizeof(uint8_t));

    for (unsigned i = 0; i < (*tree_bitmask).size(); i++)
        this->fwd_table[iface].tree_bitmask[i] = (*tree_bitmask)[i];

    this->fwd_table[iface].tree_bitmask_size = (*tree_bitmask).size();    

    // next hop information: remote router pointer and iface
    this->fwd_table[iface].next_hop.router = next_hop_router;
    this->fwd_table[iface].next_hop.iface = next_hop_iface;

    return 0;
}

RID_Router::nw_address RID_Router::get_next_hop(uint8_t iface) {
    return this->fwd_table[iface].next_hop;
}

RID_Router::fwd_table_row * RID_Router::get_fwd_table() {
    return this->fwd_table;
}

uint64_t RID_Router::get_fwd_table_size() { 
    return this->fwd_table_size; 
}

int RID_Router::get_num_entries(uint8_t iface) {
    return (int) this->fwd_table[iface].num_entries;
}

uint8_t RID_Router::get_iface_num() { 
    return this->iface_num; 
}

uint8_t RID_Router::get_f_max() { 
    return this->f_max; 
}

int RID_Router::get_tree_bitmask_size(uint8_t iface) {
    return this->fwd_table[iface].tree_bitmask_size;
}

uint8_t * RID_Router::get_tree_bitmask(uint8_t iface) {
    return this->fwd_table[iface].tree_bitmask;
}

bool RID_Router::is_iface_on_ptree(int iface, uint8_t * tree_bitmask) {

    std::cout << "RID_Router::is_iface_on_ptree() : [INFO] is " << iface 
        << " on prefix tree?";

    for (int i = 0; i < this->fwd_table[iface].tree_bitmask_size; i++) {

        uint8_t iface_tree_bitmask = this->fwd_table[iface].tree_bitmask[i];

        std::cout << "\n\t tree_bitmask vs. iface_tree_bitmask : " 
            << (int) tree_bitmask[i] << " <-> " << (int) iface_tree_bitmask;

        if (tree_bitmask[i] & iface_tree_bitmask) {
            std::cout << "\n\t IT IS!" << std::endl;
            return true;
        }
    }

    std::cout << "\n\t IT'S NOT..." << std::endl;
    return false;
}

int RID_Router::calc_ptree_iface_probs() {

    this->ptree_iface_probs.clear();

    for (int f = 0; f < this->f_max; f++) {

        this->ptree_iface_probs.push_back(0.0);

        for (int iface = IFACE_LOCAL; iface < this->iface_num; iface++) {
            // only gather values if the iface is possibly on ptree
            if (!(this->iface_on_ptree[iface]))
                continue;

            this->ptree_iface_probs[f] += this->fwd_table[iface].iface_proportion * this->fwd_table[iface].f_distribution[f];
        }
    }

    return 0;
}

int RID_Router::calc_egress_ptree_probs(
    uint8_t mode,
    uint8_t iface,
    __float080 * log_prob_fp_not_larger_than,
    __float080 * log_prob_not_fp) {

    if (mode == EGRESS_PTREE_PROB_MODE_GLOBAL) {

        for (uint8_t f = 0; f <= this->f_max; f++)
            this->egress_ptree_prob[iface][f] *= this->egress_iface_prob[iface][f];

        return 0;
    }

    // intialize egress_ptree_prob[iface][ptree_size] w/ 0.0s
    // FIXME : for now let's keep egress_ptree_prob[iface][0] = 0.0
    for (uint8_t ptree_size = 0; ptree_size <= this->f_max; ptree_size++)
        this->egress_ptree_prob[iface][ptree_size] = 0.0;

    // there are 2 events which can make iface to belong to a FP prefix tree 
    // of size p:
    //
    //  1) iface is associated w/ the a previous prefix tree of size p, i.e. 
    //     the prefix tree chosen at a previous router. this can only happen 
    //     if iface is eligible to be a continuation of the prefix tree.
    //
    //  2) iface is not associated w/ a previous prefix tree of size p BUT it 
    //     experiences a FP match of size p. in this case, we say the 
    //     request falls into a 'new' prefix tree, which starts at this router.
    //
    // as such, egress_ptree_prob[iface][p] = P(event 1) + P(event 2). also, 
    // note that P(event 2) = (1.0 - P(event 1)) * P('iface has FP size p'). 
    int iface_bitmask = this->fwd_table[iface].f_bitmask;
    for (uint8_t ptree_size = ffs(iface_bitmask); ptree_size != 0; ptree_size = ffs(iface_bitmask)) {
        iface_bitmask = iface_bitmask & (iface_bitmask - 1);

        // contribution of event 1

        // P(event 1) is the likelihood of having iface associated with an previous 
        // prefix tree. this happens if the iface is eligible. so let's calculate 
        // the probability of event 1.
        __float080 prob_iface_in_ptree = 0.0;
        if (this->iface_on_ptree[iface]) {

            // accounts for sub-event B
            prob_iface_in_ptree = this->ingress_ptree_prob[ptree_size];
            // probability of ptree_iface being in the prefix tree, i.e. sub-event C
            prob_iface_in_ptree *= (fwd_table[iface].iface_proportion * fwd_table[iface].f_distribution[ptree_size - 1]);
            prob_iface_in_ptree /= (this->ptree_iface_probs[ptree_size - 1]);

            // add the contribution of event 1 to egress_ptree_prob[iface][ptree_size]
            this->egress_ptree_prob[iface][ptree_size] += prob_iface_in_ptree;
        }

        // add contribution of P(event 2), multiplied by (1.0 - P(event 1))
        this->egress_ptree_prob[iface][ptree_size] +=
            (1.0 - prob_iface_in_ptree)
            * exp(log_prob_fp_not_larger_than[ptree_size])              // not having a FP > f
            * (1.0 - exp(log_prob_not_fp[ptree_size - 1]))              // generating a FP of size f
            * ((this->fwd_table[iface].num_entries > 0) ? 1.0 : 0.0)    // making sure that there are entries in this iface
            * ((this->fwd_table[iface].f_distribution[ptree_size - 1] > 0.0) ? 1.0 : 0.0);
    }

    return 0;
}

int RID_Router::forward(
    uint8_t request_size,
    uint8_t ingress_iface,
    int * tp_sizes,
    uint8_t * tree_bitmask,
    __float080 ingress_prob,
    __float080 * ingress_ptree_prob,
    __float080 * f_r_distribution) {

    std::cout << "RID_Router::forward() : [INFO] forwarding " << 
        "\n\tfrom : " << this->get_next_hop(ingress_iface).router->get_id() << 
        "\n\tinto : " << this->get_id() << "[" << (int) ingress_iface << "]" << std::endl;

    this->tp_sizes = tp_sizes;

    // check which ifaces are eligible to be on prefix trees
    this->iface_on_ptree.clear();
    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
        this->iface_on_ptree.push_back(is_iface_on_ptree(i, tree_bitmask));

    // pre-calculate ptree iface probabilities
    calc_ptree_iface_probs();

    // basic loop prevention : we keep a list of ifaces over which the request
    // should not be forwarded
    this->no_forwarding.clear();
    // don't foward over ingress iface
    this->no_forwarding.insert(ingress_iface);
    // don't forward over ifaces w/ no entries associated w/ them
    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
        if (this->fwd_table[i].num_entries == 0)
            this->no_forwarding.insert(i);

    std::cout << "RID_Router::forward() : [INFO] not forwarding over ifaces {";
    for (std::set<uint8_t>::iterator it = this->no_forwarding.begin();
        it != this->no_forwarding.end();
        ++it) {

        std::cout << (int) (*it) << ", ";        
    }
    std::cout << "}" << std::endl;

    // save the ingress probabilities for each diff. prefix tree size (i.e. 
    // the probability that a packet is already bound to a tree of size |F|)

    // ingress probabilities, what are these?
    //  -# ingress_prob : probability of having the request forwarded 
    //                    INTO this router (by the previous router), regardless 
    //                    of the cause. the causes can be 3:
    //                      -# true positive (in which case ingress_prob = 1.0)
    //                      -# stuck in a 'false positive' prefix tree. the 
    //                         prefixes announced in the tree can be of 
    //                         size 1, 2, 3, ..., f_max. the probability for a 
    //                         specific size is given by the ingress_ptree_prob 
    //                         array (see below).
    //                      -# fallback address
    //
    //  -# ingress_ptree_prob : probability of having the request 'stuck' in
    //                          a 'false positive' prefix tree of size ptree_size,
    //                          with 1 <= ptree_size <= f_max.
    this->ingress_prob = ingress_prob;
    std::cout << "RID_Router::forward() : P(INGRESS) = " 
        << this->ingress_prob << std::endl;

    for (uint8_t ptree_size = 0; ptree_size <= this->f_max; ptree_size++) {
        this->ingress_ptree_prob[ptree_size] = ingress_ptree_prob[ptree_size];

        std::cout << "RID_Router::forward() : P(INGRESS_PTREE[" 
            << (int) ptree_size << "]) = " << this->ingress_ptree_prob[ptree_size] << std::endl;
    }
    
    // 1) calculate distributions for the L_{i,ptree_size} random variables (RVs), 
    // for each iface i and ptree_size:
    //  -# L_{i,ptree_size} expresses the the prob of having a match of 
    //     size L (0 <= L <= f_max), at iface i, considering that the incoming
    //     request is bound to a 
    //     prefix tree of size ptree_size (0 <= ptree_size <= f_max). 
    //  -# ptree_size = 0 accounts for 'new' false positive matches, i.e. those 
    //     which happen
    //     at the current router for the 1st time in the request's lifetime.
    //  -# ptree_size > 0 accounts for false positive matches happening due to 
    //     the request being stuck to a prefix tree of size ptree_size.
    if (calc_largest_match_distributions(
        request_size, 
        tp_sizes,
        f_r_distribution) < 0) {

        return -1;
    }

    // 1.1) print L_{i,ptree_size} for all (i,ptree_size) pairs
    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++)
        this->print_lpm_pmf(iface);

    // 1.2) print the egress probabilities, egress_ptree_prob[i][f], i.e. the 
    // probability of having the request bound to a fp prefix tree of size f, 
    // over iface i
    std::cout << "RID_Router::forward() : egress ptree size probs:" << std::endl;
    for (uint8_t i = 0; i < this->iface_num; i++) {
        for (uint8_t f = 0; f <= this->f_max; f++) {

            std::cout << "\tP(EGRESS_PTREE[" << (int) i << "][" << (int) f << "]) = " 
                << this->egress_ptree_prob[i][f] << std::endl;
        }
    }

    // 2) calculate the joint probability distribution of the L RVs, i.e. 
    // L_{0,ptree_size} x L_{1,ptree_size} x ... x L_{f_max,ptree_size} 
    // this determines the likelihood of several L{i,ptree_size} combinations,
    // which in allows us to answer several questions:
    //
    //  -# how likely is it to forward the request over iface i? e.g. the 
    //     joint event {(L_{0} = 0), (L_{1} = 2), (L_{2} = 1), (L_{3} = 1)}
    //     chooses iface 1 as the largest matching iface, and its probability 
    //     adds to the probability of forwarding 
    //  -# how likely is it to forward the request over multiple ifaces? e.g. 
    //     the probability of having the event 
    //     {0, 2, 2, 1} tells us the 
    //     probability of having ifaces 1 and 2 as the largest matching 
    //     ifaces, and thus being chosen to forward the request.
    //  -# how likely is it to have no matches at all? that would be given 
    //     by the event {0, 0, 0, 0}.

    // 2.1) initializations

    // 2.1.1) initialize a single array to compute the joint prob distribution 
    // for all ptree_size
    std::map<std::string, __float080> joint_lpm_matrix;
    // __float080 * joint_lpm_matrix = NULL;
    // init_joint_lpm_pmf(&(joint_lpm_matrix));

    // 2.1.2) initialize iface_events_pmf
    for (uint8_t i = 0; i < EVENT_NUM; i++)
        this->iface_events_pmf[i] = 0.0;
    // 2.1.3) initialize egress iface probs
    for (uint8_t i = 0; i < this->iface_num; i++)
        for (uint8_t f = 0; f < (this->f_max + 1); f++)
            this->egress_iface_prob[i][f] = 0.0;

    // 2.2) compute the joint distr. considering the request is bound to 
    // a prefix tree of size ptree_size. under such conditions, the events 
    // are conditioned to the following:
    //  A) prefix tree of size ptree_size continues on an 'old' iface i
    //  B) new prefix tree, of size >= ptree_size, starts on a 'new' iface j
    for (uint8_t ptree_size = 0; ptree_size <= this->f_max; ptree_size++) {

        // 2.2.1) if the probability of having a request bound to a ptree of 
        // size ptree_size is 0.0, then - by definition - it is impossible for 
        // any next iface to be 'bound' to such a tree. thus, we skip this 
        // calculation (saves time).
        if (ingress_ptree_prob[ptree_size] == 0.0) {
            std::cout << "RID_Router::forward() : skipping prefix tree of size : " 
                << (int) ptree_size << std::endl;
            continue;
        }

        // 2.2.2) re-use the joint probability matrix for the computation of 
        // joint prob of each prefix tree size
        joint_lpm_matrix.clear();
        // clear_joint_lpm_pmf(&(joint_lpm_matrix));

        // 2.2.3) assuming the prefix tree of size ptree_size continues on 
        // iface ptree_iface...
        for (uint8_t ptree_iface = 0; ptree_iface < this->iface_num; ptree_iface++) {

            // skip ifaces on the 'no forwarding' list
            std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) ptree_iface);
            if (it != this->no_forwarding.end())
                continue;

            calc_joint_largest_match_distributions(ptree_size, ptree_iface, &joint_lpm_matrix);
        }

        // 2.3) calc the prob distribution for ifaces, for each ptree_size
        if (calc_iface_events_distributions(ptree_size, &joint_lpm_matrix) < 0)
            return -1;
    }

    std::cout << "RID_Router::forward() :"
        << "\n\t P('NLM') = " << this->iface_events_pmf[EVENT_NLM]
        << "\n\t P('LLM') = " << this->iface_events_pmf[EVENT_LLM]
        << "\n\t P('SLM') = " << this->iface_events_pmf[EVENT_SLM]
        << "\n\t P('MLM') = " << this->iface_events_pmf[EVENT_MLM]
        << std::endl;

    this->print_iface_events_pmf();
    this->print_egress_iface_prob();

    // // clean up after yourself...
    // free(joint_lpm_matrix);
    // ... aaaaand we're out of here!
    return 0;
}

__float080 fp_rate(__float080 m, __float080 n, __float080 k, __float080 c) {

    // printf("fp_rate() :"\
    //     "\n\tm = %-.5LE"\
    //     "\n\tn = %-.5LE"\
    //     "\n\tk = %-.5LE"\
    //     "\n\t|F\\R| = %-.5LE"\
    //     "\n\tP(fp) = %-.5LE\n",
    //     m, n, k, c, (__float080) pow((1.0 - exp(-((n / m) * k))), k * c));
    return pow((1.0 - exp(-((n / m) * k))), k * c);
}

/*
 * \brief   calculates the probabilities of *NOT* having *ANY* FP match of
 *          size f (f ranging from 0 to f_max) on some iface. result is given in 
 *          log form and saved on array log_prob_not_fp (passed to the function).
 *
 * \param   m                   size of Bloom filter (BF), in bit
 * \param   n                   nr. of elements encoded in BF
 * \param   k                   nr. of hash functions used in BF
 * \param   f_entries           nr. of entries associated w/ iface
 * \param   f_distribution      distribution of entry sizes associated w/ iface
 * \param   f_r_distribution    |F\R| distributions (see paper for details)
 * \param   log_prob_not_fp     array on which results are saved
 */
int RID_Router::calc_log_prob_not_fp(
    __float080 m, 
    __float080 n, 
    __float080 k,
    __float080 f_entries,               // nr. of entries for iface
    __float080 * f_distribution, 
    __float080 * f_r_distribution,
    __float080 * log_prob_not_fp) {        // function fills this array

    __float080 _n_entries = 0.0;
    __float080 _f_entries = 0.0; 
    __float080 _subtotal_f_r = 0.0; 

    // as a precaution, reset log_prob_not_fp to 0.0
    for (int f = 0; f < this->f_max; f++)
        log_prob_not_fp[f] = 0.0;

    for (int _f = 0; _f < this->f_max; _f++) {

        _f_entries = f_entries * f_distribution[_f];
        // std::cout << "RID_Router::calc_log_prob_not_fp() : _f_entries[" << _f << "] = " << _f_entries << " (out of " << f_entries << ")" << std::endl;

        if (_f == 0) {

            log_prob_not_fp[_f] += _f_entries * log(1.0 - fp_rate(m, n, k, (__float080) (_f + 1)));  

        } else {

            _subtotal_f_r = 0.0;

            for (int _jf = 0; _jf <= _f; _jf++) {
                _subtotal_f_r += f_r_distribution[_jf];
            }

            for (int _jf = 0; _jf <= _f; _jf++) {

                _n_entries = (__float080) _f_entries * (f_r_distribution[_jf] / _subtotal_f_r);
                log_prob_not_fp[_f] += _n_entries * log(1.0 - fp_rate(m, n, k, (__float080) (_jf + 1)));

                // std::cout << "RID_Router::calc_log_prob_not_fp() : log_prob_not_fp[" 
                //     << _f << "][" << _jf << "] = " 
                //     << exp(_n_entries * log(1.0 - fp_rate(m, n, k, (__float080) (_jf + 1)))) << std::endl;

                // std::cout << "RID_Router::calc_log_prob_not_fp() : log_prob_not_fp[" 
                //     << _f << "][" << _jf << "] = " 
                //     << exp(10 * log(1.0 - fp_rate(m, n, k, (__float080) (_jf + 1)))) << std::endl;
            }
        }

        // std::cout << "RID_Router::calc_log_prob_not_fp() : log_prob_not_fp[" 
        //     << _f << "] = " << exp(log_prob_not_fp[_f]) << std::endl;
    }

    return 0;
}

void RID_Router::init_lpm_pmf(uint8_t iface) {

    // // lpm_pmf is a 2D array : it may not be initialized when 
    // // init_lpm_pmf() is first called
    // if (this->lpm_pmf == NULL) {
    //     this->lpm_pmf = 
    //         (RID_Router::lpm_pmf_row **) calloc(this->iface_num, sizeof(RID_Router::lpm_pmf_row));
    // }

    // there are |R|_{max} + 1 prefix trees, starting from ptree_size = 0, 
    // meaning 'not in any prefix tree'
    this->lpm_pmf[iface] = (RID_Router::lpm_pmf_row *) calloc(this->f_max + 1, sizeof(RID_Router::lpm_pmf_row));

    for (int _p = 0; _p < this->f_max + 1; _p++) {

        // this->lpm_pmf[iface][_p] = iface;
        // this->lpm_pmf[iface][_p] = ptree_size = ptree_size;
        this->lpm_pmf[iface][_p].lpm_pmf_prob = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
        // this->lpm_pmf[iface][_p].lpm_not_in_ptree_pmf = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
    }
}

__float080 * RID_Router::get_lpm_prob(
    uint8_t iface, 
    uint8_t ptree_size) { 

    return this->lpm_pmf[iface][ptree_size].lpm_pmf_prob;
}

__float080 RID_Router::get_lpm_prob(
    uint8_t iface, 
    uint8_t ptree_size,
    uint8_t f) { 

    return this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f];
}

/*
 * \brief   prints the distribution of L_{i,p} random variables, for a single 
 *          iface
 *
 * prints a table in the following format:
 *
 * ---------------------------------------------------------------------------------
 * [|L|]                : | 0              | 1 | ... | |R|_{max} | SUM(P(L_{i,p})) |
 * ---------------------------------------------------------------------------------
 * [P(L_{i,0})]         : | P(L_{i,0} = 0) | . | ... |    ...    |      ...        |
 *     ...              : |      ...       | . | ... |    ...    |      ...        |
 * [P(L_{i,|R|_{max}})] : |      ...       | . | ... |    ...    |      ...        |      
 * ---------------------------------------------------------------------------------
 * 
 * \param   iface   iface for which the L_{i,p} table will be printed
 *
 */
void RID_Router::print_lpm_pmf(uint8_t iface) {

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->f_max + 2); i++)
        printf("-------------");

    printf("\n[|L|]         : |");

    for (uint8_t i = 0; i < (this->f_max + 1); i++)
        printf(" %-11d|", i);

    printf(" SUM P(L_{i,p}) |");

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->f_max + 2); i++)
        printf("-------------");

    __float080 _cumulative_prob = 0.0;
    __float080 _prob = 0.0;

    for (uint8_t ptree_size = 0; ptree_size < (this->f_max + 1); ptree_size++) {

        printf("\n[P(L_{%d,%d}) ] : |", iface, ptree_size);

        for (uint8_t f = 0; f < (this->f_max + 1); f++) {

            _prob = this->get_lpm_prob(iface, ptree_size, f);
            printf(" %-.5LE|", _prob);

            _cumulative_prob += _prob;
        }

        printf("     %-.5LE|", _cumulative_prob);
        _cumulative_prob = 0.0;
    }

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->f_max + 2); i++)
        printf("-------------");

    printf("\n");
}

// /*
//  * \brief   allocates memory for an array capable of holding a joint distribution 
//  *          matrix for all L_{iface, ptree_size} RVs (i.e. for all ifaces in 
//  *          the router)
//  * 
//  * \param   joint_prob_matrix   the array to initialize
//  */
// void RID_Router::init_joint_lpm_pmf(
//     __float080 ** joint_prob_matrix) {

//     *joint_prob_matrix = (__float080 *) calloc(pow((this->f_max + 1), this->iface_num), sizeof(__float080));
// }

// /*
//  * \brief   sets the joint prob distribution array to 0.0s. this is just a 
//  *          wrapper for memset().
//  * 
//  * \param   joint_prob_matrix   the array to clear
//  */
// void RID_Router::clear_joint_lpm_pmf(__float080 ** joint_prob_matrix) {


//     memset(*joint_prob_matrix, 0, pow((this->f_max + 1), this->iface_num) * sizeof(__float080));
// }

/*
 * \brief   adds (as in 'sum', '+') a prob. value to the joint_prob_matrix 
 *          position indexed by iface_pivots. in other words this function does 
 *          joint_prob_matrix[iface_pivots] += value
 *
 * multi-dimensional arrays in C/C++ are stored in row major order. what does 
 * it mean? e.g. say we have a 3D matrix, 'the_matrix' of dimension D x D x D. 
 * to access element the_matrix[1][3][5], we de-reference a pointer as such: 
 * *(the_matrix + (1 * ((D)^2)) + (3 * ((D)^1)) + (5 * ((D)^0)))
 *
 * this is why we need a function to do something as simple as 
 * joint_prob_matrix[iface_pivots] += value
 * 
 * \param   joint_prob_matrix       the function adds (as in 'sum', '+') value 
 *                                  to the iface_pivots position in this array.
 * \param   iface_pivots            index in joint_prob_matrix to which the 
 *                                  prob. value should be added.
 * \param   value                   the prob. value to be added to the array.
 */
void RID_Router::add_joint_lpm_prob(
    std::map<std::string, __float080> * joint_prob_matrix,
    int * iface_pivots,
    __float080 value) {

    std::ostringstream iface_pivots_str("");
    for (int i = 0; i < this->iface_num; i++)
        iface_pivots_str << iface_pivots[i];

    (*joint_prob_matrix)[iface_pivots_str.str()] += value;
}

/*
 * \brief   gets a prob. value to the joint_prob_matrix 
 *          position indexed by iface_pivots. in other words this function 
 *          returns joint_prob_matrix[iface_pivots]
 *
 * multi-dimensional arrays in C/C++ are stored in row major order. what does 
 * it mean? e.g. say we have a 3D matrix, 'the_matrix' of dimension D x D x D. 
 * to access element the_matrix[1][3][5], we de-reference a pointer as such: 
 * *(the_matrix + (1 * ((D)^2)) + (3 * ((D)^1)) + (5 * ((D)^0)))
 *
 * this is why we need a function to do something as simple as 
 * returning joint_prob_matrix[iface_pivots]
 * 
 * \param   joint_prob_matrix       the function adds (as in 'sum', '+') value 
 *                                  to the iface_pivots position in this array.
 * \param   iface_pivots            index in joint_prob_matrix to which the 
 *                                  prob. value should be added.
 * \param   value                   the prob. value to be added to the array.
 */
__float080 RID_Router::get_joint_lpm_prob(    
    std::map<std::string, __float080> * joint_prob_matrix,
    int * iface_pivots) {

    std::ostringstream iface_pivots_str("");
    for (int i = 0; i < this->iface_num; i++)
        iface_pivots_str << iface_pivots[i];

    return (*joint_prob_matrix)[iface_pivots_str.str()];
}

/*
 * \brief   compute the distribution of L_{i,p} random variables
 *
 * calculate distributions for the L_{i,p} random variables (RVs). L_{i,p} 
 * expresses the the prob of having a request match a forwarding entry of size 
 * L (0 <= L <= f_max), at iface i, considering that the incoming request is 
 * bound to a FP prefix tree of size p (0 <= p <= f_max). 
 *
 * regarding FP prefix trees of size p:
 *  -# p = 0 means the request isn't bound to a prefix tree. this can happen 
 *     right at the begining or when a request finds a 'new' iface.
 *     at the current router for the 1st time in the request's lifetime.
 *  -# p > 0 means the request comes bound to a FP prefix tree of size p. in 
 *     the current router, the request may continue over that tree OR bind 
 *     to a new tree of larger size.
 *
 * the objective is to build a table like the following, one for each iface:
 *
 *              |   size of largest match for iface i   |
 *              -----------------------------------------
 * | ptree_size |  f = 0 (no match) | 1 | 2 | 3 | 4 | 5 |
 * ------------------------------------------------------
 * |          0 |         "         | " | " | " | " | " |
 * |          1 |         "         | " | " | " | " | " |
 * |          2 |         "         | " | " | " | " | " |
 * |          3 |         "         | " | " | " | " | " |
 * |          4 |         "         | " | " | " | " | " |
 * |          5 |         "         | " | " | " | " | " |
 *
 */
int RID_Router::calc_largest_match_distributions(
    uint8_t request_size, 
    int * tp_sizes, 
    __float080 * f_r_distribution) {

    __float080 * log_prob_fp_not_larger_than = (__float080 *) calloc(request_size + 1, sizeof(__float080));
    __float080 * log_prob_not_fp = (__float080 *) calloc(this->f_max, sizeof(__float080));
    // the k used in the FP rate formulas (see paper)
    __float080 k = (log(2) * ((__float080) this->bf_size)) / ((__float080) request_size);

    // calculation of L_{i,ptree_size} for iface i starts here...
    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++) {

        // if iface isn't initialized, abort
        if (this->fwd_table[iface].iface == -1) {

            fprintf(stderr, "RID_Router::calc_largest_match_distributions() : iface %d uninitialized. aborting.\n", iface);
            return -1;
        }

        // allocate memory for L_{iface,ptree_size} matrix (only for iface).
        // for coding purposes, we abbreviate L_{iface,ptree_size} as lpm_pmf.
        this->init_lpm_pmf(iface);

        // we don't forward packets over forbidden ifaces, as such the 
        // P(L_{forbidden iface,<for all ptree_size>} = 0) = 1.0
        for (int ptree_size = this->f_max; ptree_size >= 0; ptree_size--) {
            
            std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) iface);
            if (it != this->no_forwarding.end())
                this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[0] = 1.0;
        }

        // skip iface if it's in the 'no forwarding' list
        std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) iface);
        if (it != this->no_forwarding.end())
            continue;

        // FIXME : if there are no entries for this iface, a match is impossible 
        // for this [iface][ptree_size][f] combination. also, set 
        // [iface][ptree_size][0] = 1.0 (a 'no match' event is guaranteed)
        if (this->fwd_table[iface].num_entries < 1) {

            for (int ptree_size = this->f_max; ptree_size >= 0; ptree_size--)
                this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[0] = 1.0;

            continue;
        }

        // pre-calculate prob of *NOT* having *ANY* FP of sizes 0 to f 
        // in iface
        if (calc_log_prob_not_fp(
            (__float080) (this->bf_size), (__float080) request_size, k,
            (__float080) this->fwd_table[iface].num_entries,
            this->fwd_table[iface].f_distribution,
            f_r_distribution,
            log_prob_not_fp) < 0) {

            free(log_prob_not_fp);
            fprintf(stderr, "RID_Router::calc_largest_match_distributions() : couldn't retrieve log FP rates\n");

            return -1;
        }

        // pre-calculate the probabilities of *NOT* having FPs larger than some 
        // size f. this will be useful for the L_{i,ptree_size} calculations. e.g. 
        // note that L_{iface,ptree_size} = f can only happen if a match larger 
        // than f doesn't occur.
        for (int f = (this->f_max - 1); f >= 0; f--)
            log_prob_fp_not_larger_than[f] = log_prob_fp_not_larger_than[f + 1] + log_prob_not_fp[f];

        // calculate the *LOCAL* component of the egress prefix tree probabilities
        calc_egress_ptree_probs(
            EGRESS_PTREE_PROB_MODE_LOCAL,
            iface,
            log_prob_fp_not_larger_than,
            log_prob_not_fp);

        for (int ptree_size = this->f_max; ptree_size >= 0; ptree_size--) {
            for (int f = this->f_max; f >= 0; f--) {

                // if f is smaller than the max. true positive match 
                // (tp_sizes[iface]) then f *WON'T BE* the largest match for 
                // sure : notice we're *GUARANTEED* to have a match of at least 
                // size tp_sizes[iface]. 
                if (f < tp_sizes[iface]) {

                    // // f should also be smaller than the current ptree_size, 
                    // // since f >= ptree_size would guarantee a positive 
                    // // match of *AT LEAST* ptree_size
                    // // FIXME 1: this doesn't make sense...
                    // if (f < ptree_size)
                    this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 0.0;

                    // FIXME 2: you may need to do something about the egress 
                    // probabilities here...

                } else {

                    // 2 sub-cases (deterministic):
                    //  1) TP exists and its size is f
                    //  2) no TP exists for f
                    if (f == tp_sizes[iface]) {

                        // sub-case 1) TP exists and its size is f

                        // if iface is in the ptree and f < ptree_size : 
                        // since ptree_size < tp_sizes[iface], f will never 
                        // be the largest match
                        if (f < ptree_size) {

                            // FIXME: this seems hack-ish : here i check if 
                            // iface is associated with any entries. why set 
                            // it to 1.0 if the iface has no entries? is this  
                            // even possible to reach?
                            this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 
                                (__float080) ((this->fwd_table[iface].num_entries > 0) ? 0.0 : 1.0);

                        } else {

                            // if f >= ptree_size, the largest match will *AT 
                            // LEAST* be f.
                            this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 
                                exp(log_prob_fp_not_larger_than[f]);
                        }

                    } else {

                        // sub-case 2) no TP exists for f
                        if (f > 0) {

                            // // if f is less than the prefix tree size, then 
                            // // any match will be larger than f
                            // if (f < ptree_size) {

                            //     this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 0.0;

                            // } else if (f == ptree_size) {
                            if (f == ptree_size) {

                                // the prob of this event can be calculated as:
                                //  * not having fp matches larger than f : log_prob_fp_not_larger_than[f]
                                //  * AND having a fp match of size f : 1.0 - exp(log_prob_not_fp[f - 1])
                                this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 
                                    exp(log_prob_fp_not_larger_than[f]) 
                                    * (1.0 - exp(log_prob_not_fp[f - 1])) 
                                    * ((this->fwd_table[iface].num_entries > 0) ? 1.0 : 0.0);

                            } else {

                                // if the iface is potentially on ptree, then it can't have 
                                // a match larger than the prefix tree. if 
                                // the iface is def. not on ptree, it can be the start of 
                                // a new prefix tree, and as such fp matches 
                                // larger than f are possible.
                                if (this->iface_on_ptree[iface] && ptree_size > 0) {

                                    this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 0.0;
                                
                                } else {

                                    // the prob of this event can be calculated as:
                                    //  * not having fps larger than f : log_prob_fp_not_larger_than[f]
                                    //  * AND having a fp of size f : 1.0 - exp(log_prob_not_fp[f - 1])
                                    this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 
                                        exp(log_prob_fp_not_larger_than[f]) 
                                        * (1.0 - exp(log_prob_not_fp[f - 1])) 
                                        * ((this->fwd_table[iface].num_entries > 0) ? 1.0 : 0.0);
                                }
                            }

                        } else {

                            // FIXME : this also seems like something we could 
                            // cut
                            if (f < ptree_size) {

                                this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 
                                    (__float080) ((this->fwd_table[iface].num_entries > 0) ? 0.0 : 1.0);

                            } else {

                                this->lpm_pmf[iface][ptree_size].lpm_pmf_prob[f] = 
                                    exp(log_prob_fp_not_larger_than[f]);
                            }
                        }
                    }
                }
            }
        }
    }

    free(log_prob_fp_not_larger_than);
    free(log_prob_not_fp);

    return 0;
}

__float080 RID_Router::calc_log_joint_largest_match_prob(
    uint8_t ptree_size,
    uint8_t ptree_iface, 
    int * iface_pivots) {

    __float080 log_joint_largest_match_prob = 0.0;
    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++) {

        if ((iface == ptree_iface))
            log_joint_largest_match_prob += log(this->get_lpm_prob(iface, ptree_size, iface_pivots[iface]));
        else
            log_joint_largest_match_prob += log(this->get_lpm_prob(iface, 0, iface_pivots[iface]));
    }

    return log_joint_largest_match_prob;
}

__float080 RID_Router::calc_joint_largest_match_prob_sum(
    uint8_t ptree_size,
    uint8_t ptree_iface) {

    // we represent joint events as a pivot array of size iface_num. e.g.
    //
    //        iface index : 0  1  2  3
    //                      v  v  v  v 
    //      iface_pivots = {0, 1, 2, 1} 
    //
    // represents the joint event in which iface 0 has no match (size 0), iface 
    // 1 has a match of size 1, iface 2 a match of size 2 and iface 3 a match 
    // of size 1. 
    // the idea is to cycle through all 
    // combinations, and calculate the probability of that event by multiplying 
    // the individual L_{i,ptree_size} and saving it in joint_prob_matrix.
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    // max. possible fwd entry sizes for each iface
    int * iface_pivots_max = (int *) calloc(this->iface_num, sizeof(int));
    // bitmasks saving unvisited fwd entry sizes for each iface
    int * iface_bitmasks = (int *) calloc(this->iface_num, sizeof(int));

    // initialize the iface_pivots_max and iface_bitmasks arrays
    for (int i = 0; i < this->iface_num; i++) {

        iface_bitmasks[i] = this->fwd_table[i].f_bitmask;
        int f_bitmask = this->fwd_table[i].f_bitmask;
        int max_pivot = 0;

        if ((ptree_size > 0) && (i == ptree_iface)) {
            iface_pivots_max[i] = ptree_size;
        } else {

            do {
                max_pivot = ffs(f_bitmask);
                f_bitmask = f_bitmask & (f_bitmask - 1);
            } while(ffs(f_bitmask) != 0);

            iface_pivots_max[i] = max_pivot;
        }
    }

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;
    // placeholder for the individual probabilities of joint events
    __float080 joint_prob_sum = 0.0;

    while (curr_i > -1) {

        // the inner loop does all the work
        if (curr_i == (this->iface_num - 1)) {

            iface_bitmasks[curr_i] = this->fwd_table[curr_i].f_bitmask;
            int f = 0;
            do {

                // update the f size pivot for iface curr_i
                iface_pivots[curr_i] = f;

                joint_prob_sum += exp(this->calc_log_joint_largest_match_prob(ptree_size, ptree_iface, iface_pivots));

                std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) curr_i);
                if (it != this->no_forwarding.end()) {
                    iface_pivots[curr_i] = iface_pivots_max[curr_i];
                    break;
                }

                f = ffs(iface_bitmasks[curr_i]);
                iface_bitmasks[curr_i] = iface_bitmasks[curr_i] & (iface_bitmasks[curr_i] - 1);

            } while((f != 0) && (f <= iface_pivots_max[curr_i]));
        }

        std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) curr_i);
        if ((it != this->no_forwarding.end()) && (iface_pivots[curr_i] > -1))
            iface_pivots[curr_i] = iface_pivots_max[curr_i];

        // extract the current pivot value for the curr_i level
        curr_f = iface_pivots[curr_i];

        if (++curr_f < (iface_pivots_max[curr_i] + 1)) {

            // we're at some level curr_i, and we need to increment 
            // it's pivot counter (to curr_f)
            if (curr_f == 0) {
                iface_pivots[curr_i] = curr_f;
            } else {
                iface_pivots[curr_i] = ffs(iface_bitmasks[curr_i]);
                iface_bitmasks[curr_i] = iface_bitmasks[curr_i] & (iface_bitmasks[curr_i] - 1);
            }

            // std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] (-*) [" << curr_i << "],[" << curr_f << "]" << std::endl;

            // after updating iface_pivots[curr_i], go back down 
            // to the level below.
            curr_i++;

        } else {

            // we've completed the iterations for the curr_i level. we 
            // now: 
            //  * reset iface_pivots[curr_i] to 0
            curr_f = -1;
            iface_pivots[curr_i] = curr_f;
            iface_bitmasks[curr_i] = this->fwd_table[curr_i].f_bitmask;
            //  * go up one level (i.e. to curr_i - 1) to increment 
            //      iface_pivots[curr_i - 1]
            // std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] (*-) [" << curr_i << "],[" << curr_f << "]" << std::endl;
            curr_i--;
        }
    }

    // delete the memory allocated by the function
    free(iface_pivots);
    free(iface_pivots_max);
    free(iface_bitmasks);

    std::cout << "RID_Router::calc_joint_largest_match_prob_sum() : [INFO] JOINT_PROB_SUM["
        << (int) ptree_iface << "][" << (int) ptree_size <<"] = " << joint_prob_sum << std::endl;

    return joint_prob_sum;
}

/*
 * \brief   compute the joint probability distribution of L_{i,p} x 
 *          L_{j,p} x ... x L_{n,p} RVs, between all ifaces in the router and 
 *          for a 'prefix tree' size p.
 * 
 * \param   ptree_size              the joint distribution should consider this 
 *                                  'prefix tree' size *ONLY*.
 * \param   ptree_iface             use the L_{i,ptree_size} values for ptree_iface, 
 *                                  L_{i,0} values for all other ifaces.
 * \param   joint_prob_matrix       save the joint distribution computations in 
 *                                  this array.
 *
 */
int RID_Router::calc_joint_largest_match_distributions(
    uint8_t ptree_size,
    uint8_t ptree_iface,
    std::map<std::string, __float080> * joint_prob_matrix) {

    // we represent joint events as a pivot array of size iface_num. e.g.
    //
    //        iface index : 0  1  2  3
    //                      v  v  v  v 
    //      iface_pivots = {0, 1, 2, 1} 
    //
    // represents the joint event in which iface 0 has no match (size 0), iface 
    // 1 has a match of size 1, iface 2 a match of size 2 and iface 3 a match 
    // of size 1. 
    // the idea is to cycle through all 
    // combinations, and calculate the probability of that event by multiplying 
    // the individual L_{i,ptree_size} and saving it in joint_prob_matrix.
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    // max. possible fwd entry sizes for each iface
    int * iface_pivots_max = (int *) calloc(this->iface_num, sizeof(int));
    // bitmasks saving unvisited fwd entry sizes for each iface
    int * iface_bitmasks = (int *) calloc(this->iface_num, sizeof(int));

    // initialize the iface_pivots_max and iface_bitmasks arrays
    for (int i = 0; i < this->iface_num; i++) {

        iface_bitmasks[i] = this->fwd_table[i].f_bitmask;
        int f_bitmask = this->fwd_table[i].f_bitmask;
        int max_pivot = 0;

        if ((ptree_size > 0) && (i == ptree_iface)) {
            iface_pivots_max[i] = ptree_size;
        } else {

            do {
                max_pivot = ffs(f_bitmask);
                f_bitmask = f_bitmask & (f_bitmask - 1);
            } while(ffs(f_bitmask) != 0);

            iface_pivots_max[i] = max_pivot;
        }
    }

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;

    __float080 prob_iface_in_ptree = 0.0;
    if (ptree_size > 0) {

        // if the ptree_iface is def. not on ptree, then it's impossible for it to be the 
        // continuation of a prefix tree.
        if (this->iface_on_ptree[ptree_iface]) {

            std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] ("
                << (int) ptree_iface << " ," << (int) ptree_size << ") continues...";

            // accounts for sub-event B
            prob_iface_in_ptree = this->ingress_ptree_prob[ptree_size];
            std::cout << "\n\t (event B) prob_iface_in_ptree = " << prob_iface_in_ptree;
            // probability of ptree_iface being in the prefix tree, i.e. sub-event C
            prob_iface_in_ptree *= (fwd_table[ptree_iface].iface_proportion * fwd_table[ptree_iface].f_distribution[ptree_size - 1]);
            std::cout << "\n\t (event C.1) prob_iface_in_ptree = " << prob_iface_in_ptree;
            prob_iface_in_ptree /= (this->ptree_iface_probs[ptree_size - 1]);
            std::cout << "\n\t (event C.2) prob_iface_in_ptree = " << prob_iface_in_ptree << std::endl;
        }

    } else {

        __float080 nr_valid_egress_ifaces = (__float080) (this->iface_num - this->no_forwarding.size());
        if (nr_valid_egress_ifaces < 1.0) nr_valid_egress_ifaces = 1.0;

        prob_iface_in_ptree = (1.0 / nr_valid_egress_ifaces);
    }

    // if prob_iface_in_ptree is 0.0, no need for further processing 
    if (prob_iface_in_ptree == 0.0)
        return 0;

    // placeholder for the individual probabilities of joint events
    __float080 log_prob = 0.0;
    // get the sum of probs of the joint events involved in this 
    // <ptree_size, ptree_iface> pair. 
    // FIXME : get_joint_prob_sum() has the same basic structure as 
    // calc_joint_largest_match_distributions(), and its invocation here 
    // seems overkill. unfortunately, i haven't found another way to deal 
    // with this issue.
    __float080 log_joint_prob_sum = log(calc_joint_largest_match_prob_sum(ptree_size, ptree_iface));

    while (curr_i > -1) {

        // the inner loop does all the work
        if (curr_i == (this->iface_num - 1)) {

            iface_bitmasks[curr_i] = this->fwd_table[curr_i].f_bitmask;
            int f = 0;
            do {

                // update the f size pivot for iface curr_i
                iface_pivots[curr_i] = f;

                // notice we're calculating the joint probability of the joint 
                // event, which accounts for 3 things:
                //  A) match lengths as encoded on iface_pivots
                //  B) request bound to a prefix tree of size ptree_size
                //  C) prefix tree of size ptree_size continues on ptree_iface
                log_prob = this->calc_log_joint_largest_match_prob(ptree_size, ptree_iface, iface_pivots);

                // scale the event probability by the following factors:
                //  1) probability of having a request coming into the router
                //  2) probability of having ptree_iface as the continuation 
                //     of a prefix tree of size ptree_size
                //  3) sum of all individual possible joint events, in order 
                //     to allow normalization of their probabilities
                log_prob += log(this->ingress_prob);
                if (this->tp_sizes[ptree_iface] > 0 && this->tp_sizes[ptree_iface] == ptree_size) {

                    // FIXME: this sounds hack-ish... but it works!
                    // check if iface_pivots encodes an event in which the 2 
                    // sub-events happen:
                    //  1) ptree_size == tp_sizes[ptree_iface]
                    //  2) other iface_pivots[] indexes other than ptree_iface 
                    //     is 0

                    // this represents the event of having the prefix tree of 
                    // size ptree_size continuing on ptree_iface AND having 
                    // a true positive entry of size ptree_size on ptree_iface.
                    // in this case, we should not multiply log_prob by 
                    // prob_iface_in_ptree, rather by this->ingress_prob, since 
                    // 
                    int iface_pivots_sum = 0;
                    for (uint8_t k = 0; k < this->iface_num; k++)
                        iface_pivots_sum += iface_pivots[k];

                    // if that's the case, we scale log_prob by the ingress 
                    // probability. otherwise, multiply by prob_iface_in_ptree
                    if (iface_pivots_sum != this->tp_sizes[ptree_iface])
                        log_prob += log(prob_iface_in_ptree);

                } else {

                    log_prob += log(prob_iface_in_ptree);
                }

                log_prob -= log_joint_prob_sum;

                // add the calculated probability to the matrix of joint 
                // probabilities
                __float080 prob = exp(log_prob);
                this->add_joint_lpm_prob(joint_prob_matrix, iface_pivots, prob);

                // if (prob > 0.0) {

                //     std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] JOINT_PROB" 
                //         << "(" << (int) ptree_iface << ", " << (int) ptree_size << ")";
                //     for (int k = 0; k < this->iface_num; k++)
                //         std::cout << "[" << iface_pivots[k] << "]";
                //     std::cout << " = " << prob << std::endl;

                //     std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] (SUM) JOINT_PROB";
                //     for (int k = 0; k < this->iface_num; k++)
                //         std::cout << "[" << iface_pivots[k] << "]";
                //     std::cout << " = " << this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots) << std::endl;
                // }

                std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) curr_i);
                if (it != this->no_forwarding.end()) {
                    iface_pivots[curr_i] = iface_pivots_max[curr_i];
                    break;
                }

                f = ffs(iface_bitmasks[curr_i]);
                iface_bitmasks[curr_i] = iface_bitmasks[curr_i] & (iface_bitmasks[curr_i] - 1);

            } while((f != 0) && (f <= iface_pivots_max[curr_i]));
        }

        std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) curr_i);
        if ((it != this->no_forwarding.end()) && (iface_pivots[curr_i] > -1))
            iface_pivots[curr_i] = iface_pivots_max[curr_i];

        // extract the current pivot value for the curr_i level
        curr_f = iface_pivots[curr_i];

        if (++curr_f < (iface_pivots_max[curr_i] + 1)) {

            // we're at some level curr_i, and we need to increment 
            // it's pivot counter (to curr_f)
            if (curr_f == 0) {
                iface_pivots[curr_i] = curr_f;
            } else {
                iface_pivots[curr_i] = ffs(iface_bitmasks[curr_i]);
                iface_bitmasks[curr_i] = iface_bitmasks[curr_i] & (iface_bitmasks[curr_i] - 1);
            }

            // std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] (-*) [" << curr_i << "],[" << curr_f << "]" << std::endl;

            // after updating iface_pivots[curr_i], go back down 
            // to the level below.
            curr_i++;

        } else {

            // we've completed the iterations for the curr_i level. we 
            // now: 
            //  * reset iface_pivots[curr_i] to 0
            curr_f = -1;
            iface_pivots[curr_i] = curr_f;
            iface_bitmasks[curr_i] = this->fwd_table[curr_i].f_bitmask;
            //  * go up one level (i.e. to curr_i - 1) to increment 
            //      iface_pivots[curr_i - 1]
            // std::cout << "RID_Router::calc_joint_largest_match_distributions() : [INFO] (*-) [" << curr_i << "],[" << curr_f << "]" << std::endl;
            curr_i--;
        }
    }

    // delete the memory allocated by the function
    free(iface_pivots);
    free(iface_pivots_max);
    free(iface_bitmasks);

    return 0;
}

__float080 RID_Router::calc_cumulative_prob(
    uint8_t mode,
    uint8_t fixed_iface, 
    uint8_t fixed_iface_size,
    std::map<std::string, __float080> * joint_prob_matrix) {

    printf("RID_Router::calc_cumulative_prob() : [fixed_iface : %d][fixed_iface_size : %d][MODE : %s]\n", 
        fixed_iface, fixed_iface_size, CUMULATIVE_PROB_MODE_STR[mode]);
    __float080 prob = 0.0;
    __float080 cumulative_prob = 0.0;

    // depending on the way RID routers handle multiple matches, we 
    // set the probability of choosing the iface differently.
    //
    //  1) P(iface, MMH_FLOOD)  : we consider fixed_iface is chosen even 
    //                            if other ifaces share the largest 
    //                            match w/ it. e.g. event {0, 5, 5}.
    //  2) P(iface, MMH_RANDOM) : let's forget about the random...
    //  3) p(iface, -)          : only account for 'clean' single 
    //                            largest match for fixed_iface. e.g. 
    //                            if fixed_iface = 1, event {0, 5, 4},
    //                            {0, 5, 5} isn't valid. 
    bool distribute_probs = false;
    if (this->mm_mode == MMH_FLOOD && mode == MODE_EI_INCLUSIVE) {
        distribute_probs = true;
    } else if (this->mm_mode != MMH_FLOOD && mode == MODE_EI_EXCLUSIVE) {
        distribute_probs = true;
    }

    // we represent joint events as a pivot array of size iface_num. e.g.
    //
    //        iface index : 0  1  2  3
    //                      v  v  v  v 
    //      iface_pivots = {0, 1, 2, 1} 
    //
    // represents the joint event in which iface 0 has no match (size 0), iface 
    // 1 has a match of size 1, iface 2 a match of size 2 and iface 3 a match 
    // of size 1. 
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    // bitmasks saving unvisited fwd entry sizes for each iface
    int * iface_bitmasks = (int *) calloc(this->iface_num, sizeof(int));
    // max. possible fwd entry sizes for each iface
    int * iface_pivots_max = (int *) calloc(this->iface_num, sizeof(int));

    // initialize the iface_bitmasks
    for (int i = 0; i < this->iface_num; i++)
        iface_bitmasks[i] = this->fwd_table[i].f_bitmask;

    // initialize the iface_pivots_max array for f
    if (mode == MODE_EI_INCLUSIVE) {

        for (uint8_t k = 0; k < fixed_iface; k++)
            iface_pivots_max[k] = (fixed_iface_size - 1);
        for (uint8_t k = fixed_iface; k < this->iface_num; k++) {

            int max_pivot = 0;
            int f_bitmask = this->fwd_table[k].f_bitmask;

            do {
                max_pivot = ffs(f_bitmask);
                f_bitmask = f_bitmask & (f_bitmask - 1);
            } while(ffs(f_bitmask) != 0);

            iface_pivots_max[k] = std::min(max_pivot, (int) fixed_iface_size);
        }

    } else {

        for (uint8_t k = 0; k < this->iface_num; k++)
            iface_pivots_max[k] = (fixed_iface_size - 1);
    }

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;

    // indicator variable which tells direction of increments in iface_pivots 
    // values
    // FIXME: this is sort of hacky, but, well... it works ?
    uint8_t going_right = 0x00;
    // iface index in iface_pivots will *ALWAYS* be set to f
    iface_pivots[fixed_iface] = fixed_iface_size;

    while (curr_i > -1) {

        // the inner loop does all the work
        if (curr_i == (this->iface_num - 1)) {

            // use iface_bitmasks to cycle through the fwd entry sizes 
            // which actually have entries, instead of using a for() cycle 
            // which goes through all indexes.
            // FIXME: it's not clear if ffs() actually provides any benefit 
            // regarding this
            iface_bitmasks[curr_i] = this->fwd_table[curr_i].f_bitmask;
            int f = 0;
            do {

                // update the f size pivot for iface curr_i
                iface_pivots[curr_i] = f;
                // fixed_iface index in iface_pivots will *ALWAYS* be set to fixed_iface_size
                iface_pivots[fixed_iface] = fixed_iface_size;

                // if MODE_EI_INCLUSIVE, make sure we don't go over the limits 
                // of ifaces set in iface_pivots_max, which in turn avoids 
                // summing the probabilities of the same event twice
                if ((curr_i != fixed_iface) && iface_pivots[curr_i] > iface_pivots_max[curr_i]) {
                    
                    // jump off the do-while() cycle
                    break;

                } else {

                    // extract the probability of the joint event encoded in 
                    // iface pivots
                    prob = this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);

                    std::cout << "RID_Router::calc_cumulative_prob() : [INFO] CUMULATIVE_PROB" 
                        << "(" << (int) fixed_iface << ", " << (int) fixed_iface_size << ") : ";
                    for (int k = 0; k < this->iface_num; k++)
                        std::cout << "[" << iface_pivots[k] << "]";
                    std::cout << " = " << prob << std::endl;
                }

                // distribute the probabilities over the egress iface probabilities
                if (distribute_probs) {
                    for (int k = 0; k < this->iface_num; k++) {
                        if ((iface_pivots[k]) > 0 && (iface_pivots[k] == fixed_iface_size))
                            this->egress_iface_prob[k][iface_pivots[k]] += prob;
                    }
                }

                // add probability to cumulative probability
                cumulative_prob += prob;

                // extract the next valid fwd entry size from the bitmask
                f = ffs(iface_bitmasks[curr_i]);
                iface_bitmasks[curr_i] = iface_bitmasks[curr_i] & (iface_bitmasks[curr_i] - 1);

            } while((f != 0) && (iface_pivots[curr_i] <= iface_pivots_max[curr_i]) && (curr_i != fixed_iface));

            // if the current iface is fixed_iface, we update the size 
            // of the iface at the upper level (i.e. at index (iface - 1) in 
            // iface_pivots). use this to make sure we go left on the pivot 
            // array.
            if (curr_i == fixed_iface)
                going_right = 0x00;
        }

        if ((curr_i == fixed_iface) && going_right) {
            
            curr_f = 0;
            // std::cout << "(->) curr_i : " << curr_i << ", curr_f : " << curr_f << ", max : " << iface_pivots_max[curr_i] << std::endl;

        } else {

            curr_f = iface_pivots[curr_i];
            // std::cout << "(<-) curr_i : " << curr_i << ", curr_f : " << curr_f << ", max : " << iface_pivots_max[curr_i] << std::endl;
        }

        // increment the current fwd entry size for curr_i
        curr_f += 1;

        // if MODE_EI_INCLUSIVE, make sure we don't go over the limits 
        // of ifaces set in iface_pivots_max, which in turn avoids 
        // summing the probabilities of the same event twice
        if (curr_f <= iface_pivots_max[curr_i]) {

            // // FIXME: if we're only looking at EI events, we only need to consider 
            // // those for which iface_pivots[IFACE_LOCAL] = 0. as soon as we 
            // // increment iface_pivots[IFACE_LOCAL] to 1, get off the cycle.
            // if ((mode == MODE_EI_EXCLUSIVE || mode == MODE_EI_INCLUSIVE) && curr_i == IFACE_LOCAL)
            //     break;

            // we're at some iface pivot index curr_i, and we need to increment 
            // it's pivot value (to curr_f)
            if (curr_i != fixed_iface) {

                if (curr_f == 0) {
                    iface_pivots[curr_i] = curr_f;
                } else {
                    iface_pivots[curr_i] = ffs(iface_bitmasks[curr_i]);

                    if (iface_pivots[curr_i] == 0)
                        iface_pivots[curr_i] = iface_pivots_max[curr_i] + 1;

                    iface_bitmasks[curr_i] = iface_bitmasks[curr_i] & (iface_bitmasks[curr_i] - 1);
                }
            }

            // std::cout << "(-*) curr_i : " << curr_i << ", curr_f : " << curr_f << ", max : " << iface_pivots_max[curr_i] << std::endl;

            if ((curr_i != fixed_iface) && iface_pivots[curr_i] > iface_pivots_max[curr_i])
                continue;

            // after updating iface_pivots[curr_i], go back to the right, i.e. 
            // deeper in the cycle
            curr_i++;
            going_right = 0xFF;

        } else {

            // we've completed the iterations for the curr_i level. we 
            // now: 
            //  * reset iface_pivots[curr_i] to 0
            if (curr_i != fixed_iface) {

                curr_f = -1;
                iface_pivots[curr_i] = curr_f;
                iface_bitmasks[curr_i] = this->fwd_table[curr_i].f_bitmask;
            }

            // order the function go up one level (i.e. to curr_i - 1) to increment 
            // iface_pivots[curr_i - 1]. 
            // e.g. we would go from [5][0][-1] -> [5][1][-1]
            // this only happens in the next cycle though
            if (--curr_i == IFACE_LOCAL && mode == MODE_LI)
                break;

            // std::cout << "(*-) curr_i : " << curr_i << ", curr_f : " << curr_f << ", max : " << iface_pivots_max[curr_i] << std::endl;

            going_right = 0x00;
        }
    }

    free(iface_pivots);
    free(iface_bitmasks);
    free(iface_pivots_max);

    return cumulative_prob;
}

int RID_Router::calc_iface_events_distributions(
    int ptree_size,
    std::map<std::string, __float080> * joint_prob_matrix) {

    // // LOCAL MATCH (LLM event)
    // int local_iface_bitmask = this->fwd_table[IFACE_LOCAL].f_bitmask;
    // for (uint8_t fixed_iface_size = ffs(local_iface_bitmask); fixed_iface_size != 0; fixed_iface_size = ffs(local_iface_bitmask)) {

    //     local_iface_bitmask = local_iface_bitmask & (local_iface_bitmask - 1);
    //     this->iface_events_pmf[EVENT_LLM] += this->calc_cumulative_prob(
    //                                                     MODE_EI_INCLUSIVE,
    //                                                     IFACE_LOCAL,
    //                                                     fixed_iface_size,
    //                                                     joint_prob_matrix);
    // }

    // std::cout << "RID_Router::calc_iface_events_distributions() : [INFO]"
    //     << "\n\t P('EVENT_LLM' [ptree_size = " << (int) ptree_size << "]') = " << this->iface_events_pmf[EVENT_LLM]
    //     << std::endl;

    // NO MATCHES : just call get_joint_lpm_prob() with iface_pivots = {0, 0, ..., 0}
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    this->iface_events_pmf[EVENT_NLM] += this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);
    // clean up allocated memory
    free(iface_pivots);

    // MULTIPLE & SINGLE IFACE MATCHES 
    __float080 prob_single = 0.0;
    __float080 prob_multiple = 0.0;

    for (uint8_t fixed_iface = 0; fixed_iface < this->iface_num; fixed_iface++) {

        int iface_bitmask = this->fwd_table[fixed_iface].f_bitmask;
        for (uint8_t fixed_iface_size = ffs(iface_bitmask); fixed_iface_size != 0; fixed_iface_size = ffs(iface_bitmask)) {
            iface_bitmask = iface_bitmask & (iface_bitmask - 1);


            // probability of having the largest match on iface, assuming 
            // other ifaces have matches *smaller* than iface. we call this 
            // MODE_EI_EXCLUSIVE
            prob_single = this->calc_cumulative_prob(
                                    MODE_EI_EXCLUSIVE,
                                    fixed_iface,
                                    fixed_iface_size,
                                    joint_prob_matrix);

            // in some cases, we may want the probability of having the largest 
            // match on iface, even if together with any other iface (e.g. if 
            // a flooding multiple match resolution strategy is in place). we 
            // call this MODE_EI_INCLUSIVE.
            prob_multiple = this->calc_cumulative_prob(
                                    MODE_EI_INCLUSIVE,
                                    fixed_iface,
                                    fixed_iface_size,
                                    joint_prob_matrix);

            // keep track of the total EI probability as well (sanity check)
            if (prob_single > 0.0) {

                std::cout << "RID_Router::calc_iface_events_distributions() : [INFO]"
                    << "\n\t P('EGRESS_IFACE_EXCLUSIVE' [" << (int) fixed_iface << "][" << (int) fixed_iface_size << "]') = " << prob_single
                    << std::endl;
            }

            if (prob_multiple > 0.0) {

                std::cout << "RID_Router::calc_iface_events_distributions() : [INFO]"
                    << "\n\t P('EGRESS_IFACE_INCLUSIVE' [" << (int) fixed_iface << "][" << (int) fixed_iface_size << "]') = " << prob_multiple
                    << std::endl;
            }

            // save the probability of a 'clean' single link match (SLM) event
            if (fixed_iface > 0)
                this->iface_events_pmf[EVENT_SLM] += prob_single;
            else
                this->iface_events_pmf[EVENT_LLM] += prob_multiple;

            // save the probability of multiple link match (MLM) events
            this->iface_events_pmf[EVENT_MLM] += (prob_multiple - prob_single);
            std::cout << "RID_Router::calc_iface_events_distributions() : [INFO] P('EVENT_MLM' [ptree_size = " << ptree_size << "]) = " 
                << this->iface_events_pmf[EVENT_MLM] << std::endl;

            // // keep track of the total EI probability as well (sanity check)
            // if (prob_single > 0.0 || prob_multiple > 0.0) {

            //     std::cout << "RID_Router::calc_iface_events_distributions() : [INFO] "
            //         << "\n\t P('EGRESS_IFACE[" << (int) iface << "][" << (int) f << "]') = " << this->egress_iface_prob[iface][f]
            //         << std::endl;
            // }
        }
    }

    return 0;
}

__float080 * RID_Router::get_egress_ptree_prob(uint8_t iface) {

    return this->egress_ptree_prob[iface];
}

__float080 RID_Router::get_iface_events_prob(uint8_t event) {

    if (event < 0 || event > EVENT_SLM) {

        fprintf(stderr, "RID_Router::get_iface_events_prob() : unknown event nr. (%d)\n", event);

        return 0.0;
    }

    return this->iface_events_pmf[event];
}

__float080 RID_Router::get_egress_iface_prob(uint8_t iface) {

    __float080 iface_prob = 0.0;

    for (int f = 0; f < (this->f_max + 1); f++)
        iface_prob += this->egress_iface_prob[iface][f];

    return iface_prob;
}

__float080 * RID_Router::get_egress_iface_probs(uint8_t iface) {

    return this->egress_iface_prob[iface];
}

/*
 * \brief   prints the probabilities of iface events for the router
 *
 * prints a table in the following format:
 *
 * ---------------------------------------------------------
 * [EVENT e]    : | [NIS] | [MIS] | [LI] | [EI] | SUM P(e) |
 * ---------------------------------------------------------
 * [P(e)]       : |   .   |   .   |  ..  |  ..  |    ..    |
 * ---------------------------------------------------------
 *
 */
void RID_Router::print_iface_events_pmf() {

    printf("\n-----------------");
    for (uint8_t event = 0; event < EVENT_NUM + 1; event++)
        printf("------------");

    printf("\n[EVENT] :  |");

    for (int event = 0; event < EVENT_NUM; event++)
        printf(" %-11s|", EVENT_STR[event]);

    printf("  SUM P(e)  |");

    printf("\n-----------------");
    for (uint8_t event = 0; event < EVENT_NUM + 1; event++)
        printf("------------");

    __float080 prob = 0.0;
    __float080 prob_total = 0.0;

    printf("\n[P(e)] :   |");

    for (uint8_t event = 0; event < EVENT_NUM; event++) {

        prob = this->iface_events_pmf[event];
        prob_total += prob;

        printf(" %-.5LE|", prob);
    }

    printf(" %-.5LE|", prob_total);

    printf("\n-----------------");
    for (uint8_t event = 0; event < EVENT_NUM + 1; event++)
        printf("------------");

    printf("\n");
}

/*
 * \brief   prints the distribution of egress size probabilities, for all 
 *          egress ifaces in the router (iface > IFACE_LOCAL).
 *
 * prints a table in the following format:
 *
 * ---------------------------------------------------------------------
 * [|L_i|]          : | 1          | 2 | ... | |R|_{max} | SUM(P(L_i)) |
 * ---------------------------------------------------------------------
 * [P(L_1)]         : | P(L_i = 1) | . | ... |     .     |       .     |
 * [P(L_2)]         : |     ...    | . | ... |     .     |       .     |
 *     ...          : |     ...    | . | ... |     .     |       .     |
 * [P(L_n)]         : |     ...    | . | ... |     .     |       .     |      
 * ---------------------------------------------------------------------
 *
 */
void RID_Router::print_egress_iface_prob() {

    printf("\n------------");
    for (uint8_t _f = 0; _f < this->f_max + 2; _f++)
        printf("-------------");

    printf("\n[|L_i|]  : |");

    for (uint8_t _f = 0; _f < this->f_max + 1; _f++)
        printf(" %-11d|", _f);

    printf(" SUM P(L_i) |");

    printf("\n------------");
    for (uint8_t _f = 0; _f < this->f_max + 2; _f++)
        printf("-------------");

    __float080 _cumulative_prob = 0.0;
    __float080 _prob = 0.0;

    for (uint8_t _iface = IFACE_LOCAL + 1; _iface < (this->iface_num); _iface++) {

        printf("\n[P(L_%d)] : |", _iface);

        for (uint8_t _f = 0; _f < this->f_max + 1; _f++) {

            _prob = this->egress_iface_prob[_iface][_f];
            printf(" %-.5LE|", _prob);

            _cumulative_prob += _prob;
        }

        printf(" %-.5LE|", _cumulative_prob);
        _cumulative_prob = 0.0;
    }

    printf("\n------------");
    for (uint8_t _f = 0; _f < this->f_max + 2; _f++)
        printf("-------------");

    printf("\n");   
}

/*
 * \brief   prints the distribution of ptree size probabilities, for all 
 *          egress ifaces in the router (iface > IFACE_LOCAL).
 *
 * prints a table in the following format:
 *
 * ---------------------------------------------------------------------
 * [|L_i|]          : | 1          | 2 | ... | |R|_{max} | SUM(P(L_i)) |
 * ---------------------------------------------------------------------
 * [P(L_1)]         : | P(L_i = 1) | . | ... |     .     |       .     |
 * [P(L_2)]         : |     ...    | . | ... |     .     |       .     |
 *     ...          : |     ...    | . | ... |     .     |       .     |
 * [P(L_n)]         : |     ...    | . | ... |     .     |       .     |      
 * ---------------------------------------------------------------------
 *
 */
void RID_Router::print_egress_ptree_prob() {

    printf("\n------------");
    for (uint8_t _f = 0; _f < this->f_max + 2; _f++)
        printf("-------------");

    printf("\n[|L_i|]  : |");

    for (uint8_t _f = 0; _f < this->f_max + 1; _f++)
        printf(" %-11d|", _f);

    printf(" SUM P(L_i) |");

    printf("\n------------");
    for (uint8_t _f = 0; _f < this->f_max + 2; _f++)
        printf("-------------");

    __float080 _cumulative_prob = 0.0;
    __float080 _prob = 0.0;

    for (uint8_t _iface = IFACE_LOCAL + 1; _iface < (this->iface_num); _iface++) {

        printf("\n[P(L_%d)] : |", _iface);

        for (uint8_t _f = 0; _f < this->f_max + 1; _f++) {

            _prob = this->egress_ptree_prob[_iface][_f];
            printf(" %-.5LE|", _prob);

            _cumulative_prob += _prob;
        }

        printf(" %-.5LE|", _cumulative_prob);
        _cumulative_prob = 0.0;
    }

    printf("\n------------");
    for (uint8_t _f = 0; _f < this->f_max + 2; _f++)
        printf("-------------");

    printf("\n");   
}   