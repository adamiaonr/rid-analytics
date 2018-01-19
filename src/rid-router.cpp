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
    uint8_t f_max,
    uint16_t bf_size,
    int mm_mode) {

    if (this->initialized)
        return 0;

    this->id = router_id;
    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    this->f_max = f_max;
    this->bf_size = bf_size;

    // initialize forwarding table
    this->fwd_table = (RID_Router::fwd_table_row *) calloc(iface_num, sizeof(RID_Router::fwd_table_row)); 
    // set fwd_table[i].iface to -1, to mark iface i as not filled in the table
    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
        this->fwd_table[i].iface = -1;

    // largest match prob. (lmp) array
    this->lmp = (RID_Router::lmp_row **) calloc(this->iface_num, sizeof(RID_Router::lmp_row));
    // avoid repeated calculation of lmps: set lmp_calculated to true when the 
    // lmp array is filled
    this->lmp_calculated = false;

    // iface event probabilities. events are: 
    //  - EVENT_NLM: no link matches
    //  - EVENT_MLM: multiple link matches
    //  - EVENT_LLM: local link match
    //  - EVENT_SLM: single link match (other than local)
    this->iface_events_prob = (__float080 *) calloc(EVENT_NUM, sizeof(__float080));

    // initialize the ingress_ptree_prob[p] array. holds the prob of having 
    // request entering a router on a prefix tree of size p.
    this->ingress_ptree_prob = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
    // initialize the egress_ptree_prob[i][p] array. holds the prob of having 
    // request leaving a router on iface i, on a prefix tree of size p.
    this->egress_ptree_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_ptree_prob[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    // initialize the egress_iface_prob[i][f] array. holds the probability 
    // of having a request leaving a router over iface i, due to a match of 
    // size f.
    this->egress_iface_prob = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_iface_prob[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    this->initialized = true;
    this->mm_mode = mm_mode;

    return 0;
}

RID_Router::RID_Router(
    std::string router_id,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max,
    uint16_t bf_size,
    int mm_mode) {

    this->init(router_id, fwd_table_size, iface_num, f_max, bf_size, mm_mode);
}

RID_Router::~RID_Router() {

    // free forwarding table
    for (uint8_t iface = 0; iface < this->iface_num; iface++)
        free(this->fwd_table[iface].f_distribution);
    free(this->fwd_table);

    // free lmp
    for (uint8_t iface = 0; iface < this->iface_num; iface++) {

        if (this->lmp[iface] == NULL)
            continue;

        for (uint8_t _f = 0; _f < this->f_max; _f++)
            free(this->lmp[iface][_f].lmp);

        free(this->lmp[iface]);
    }

    free(this->iface_events_prob);

    for (uint8_t iface = IFACE_LOCAL + 1; iface < this->iface_num; iface++)
        free(this->egress_ptree_prob[iface]);

    for (uint8_t iface = IFACE_LOCAL + 1; iface < this->iface_num; iface++)
        free(this->egress_iface_prob[iface]);
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

        // FIXME: this special case handling sounds hack-ish... but it works!
        if (this->tp_sizes[iface] > 0) {

            // int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
            // for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
            //     iface_pivots[i] = this->tp_sizes[i];
            // this->egress_ptree_prob[iface][0] += exp(this->calc_log_forwarding_event_prob(this->tp_sizes[iface], iface, iface_pivots));
            this->egress_ptree_prob[iface][0] += this->prob_true_negative;

            std::cout << "RID_Router::calc_egress_ptree_probs() : [INFO] added ptree_size 0 prob" 
                << " = " << this->egress_ptree_prob[iface][0] << std::endl;

            std::cout << "RID_Router::calc_egress_ptree_probs() : [INFO] shouldn't it be" 
                << " = " << this->prob_true_negative << " ?" << std::endl;

            // free(iface_pivots);
        }

        __float080 egress_ptree_prob_sum = 0.0;
        for (uint8_t f = 0; f <= this->f_max; f++)
            egress_ptree_prob_sum += this->egress_ptree_prob[iface][f];

        if (egress_ptree_prob_sum == 0.0)
            return 0;

        for (uint8_t f = 0; f <= this->f_max; f++)
            this->egress_ptree_prob[iface][f] /= egress_ptree_prob_sum;

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
    uint8_t * tree_bitmask,     // FIXME: shouldn't this be moved to RID_Router?
    __float080 ingress_prob,
    __float080 * ingress_ptree_prob,
    __float080 * f_r_distribution) {

    std::cout << "RID_Router::forward() : [INFO] forwarding " << 
        "\n\tfrom : " << this->get_next_hop(ingress_iface).router->get_id() << 
        "\n\tinto : " << this->get_id() << "[" << (int) ingress_iface << "]" << std::endl;

    this->tp_sizes = tp_sizes;

    // check which ifaces are eligible to be on prefix trees (of any size)
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
    // the probability that a packet is already bound to a tree of size p)

    // ingress probabilities, what are these?
    //  - ingress_prob : probability of having the request forwarded 
    //                   INTO this router by the previous router, regardless 
    //                   of the cause.
    //
    //  - ingress_ptree_prob : probability of having the request coming in while 
    //                         'stuck' in a 'false positive' prefix tree of 
    //                         size p, with 1 <= p <= f_max.
    this->ingress_prob = ingress_prob;
    std::cout << "RID_Router::forward() : P(INGRESS) = " 
        << this->ingress_prob << std::endl;

    for (uint8_t ptree_size = 0; ptree_size <= this->f_max; ptree_size++) {
        this->ingress_ptree_prob[ptree_size] = ingress_ptree_prob[ptree_size];
        std::cout << "RID_Router::forward() : P(INGRESS_PTREE[" 
            << (int) ptree_size << "]) = " << this->ingress_ptree_prob[ptree_size] << std::endl;
    }
    
    // 1) calculate distributions for the L_{i,p} random variables (RVs), 
    // for each iface i and p
    if (calc_lmp(request_size, tp_sizes, f_r_distribution) < 0)
        return -1;

    // 1.1) print L_{i,p} for all (i,p) pairs
    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++)
        this->print_lmp(iface);

    // 1.2) print egress_ptree_prob[i][p], i.e. the 
    // probability of having the request bound to a fp prefix tree of size p, 
    // over iface i
    std::cout << "RID_Router::forward() : egress ptree probs:" << std::endl;
    for (uint8_t i = 0; i < this->iface_num; i++) {
        for (uint8_t p = 0; p <= this->f_max; p++) {
            std::cout << "\tP(EGRESS_PTREE[" << (int) i << "][" << (int) p << "]) = " 
                << this->egress_ptree_prob[i][p] << std::endl;
        }
    }

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
    std::map<std::string, __float080> joint_lmp;

    // 2.1.2) initialize iface_events_prob
    for (uint8_t i = 0; i < EVENT_NUM; i++)
        this->iface_events_prob[i] = 0.0;
    // 2.1.3) initialize egress iface probs
    for (uint8_t i = 0; i < this->iface_num; i++)
        for (uint8_t f = 0; f < (this->f_max + 1); f++)
            this->egress_iface_prob[i][f] = 0.0;

    // 2.2) compute the fwd event probs given that the request is bound to 
    // a fp tree of size s. under such conditions, the possible events are 
    // restricted to:
    //  - fp tree of size s continues from the previous router, OR (nor XOR)
    //  - new fp tree, of size >= s, starts at this router
    for (uint8_t ptree_size = 0; ptree_size <= this->f_max; ptree_size++) {

        // 2.2.1) if the probability of having a request bound to a fp tree of 
        // size s is 0.0, then there's no binding. 
        // thus, we skip this calculation (saves time).
        if (ingress_ptree_prob[ptree_size] == 0.0) {
            std::cout << "RID_Router::forward() : skipping prefix tree of size : " << (int) ptree_size << std::endl;
            continue;
        }

        // 2.2.2) re-use the fwd event prob matrix for each fp tree size
        joint_lmp.clear();

        // 2.2.3) fp tree of size s continues on iface i
        for (uint8_t ptree_iface = 0; ptree_iface < this->iface_num; ptree_iface++) {

            // skip ifaces on the 'no forwarding' list
            std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) ptree_iface);
            if (it != this->no_forwarding.end())
                continue;

            calc_forwarding_event_probs(ptree_size, ptree_iface, &joint_lmp);
        }

        // 2.3) calc egress iface (and events) probs (assuming fp tree size s)
        if (calc_iface_probs(ptree_size, &joint_lmp) < 0)
            return -1;
    }

    std::cout << "RID_Router::forward() :"
        << "\n\t P('NLM') = " << this->iface_events_prob[EVENT_NLM]
        << "\n\t P('LLM') = " << this->iface_events_prob[EVENT_LLM]
        << "\n\t P('SLM') = " << this->iface_events_prob[EVENT_SLM]
        << "\n\t P('MLM') = " << this->iface_events_prob[EVENT_MLM]
        << std::endl;

    for (uint8_t i = 0; i < this->iface_num; i++)
        calc_egress_ptree_probs(EGRESS_PTREE_PROB_MODE_GLOBAL, i, NULL, NULL);

    this->print_iface_events_prob();
    this->print_egress_iface_prob();

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

void RID_Router::init_lmp(uint8_t iface) {

    // // lmp is a 2D array : it may not be initialized when 
    // // init_lmp() is first called
    // if (this->lmp == NULL) {
    //     this->lmp = 
    //         (RID_Router::lmp_row **) calloc(this->iface_num, sizeof(RID_Router::lmp_row));
    // }

    // there are |R|_{max} + 1 prefix trees, starting from ptree_size = 0, 
    // meaning 'not in any prefix tree'
    this->lmp[iface] = (RID_Router::lmp_row *) calloc(this->f_max + 1, sizeof(RID_Router::lmp_row));

    for (int _p = 0; _p < this->f_max + 1; _p++) {

        // this->lmp[iface][_p] = iface;
        // this->lmp[iface][_p] = ptree_size = ptree_size;
        this->lmp[iface][_p].lmp = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
        // this->lmp[iface][_p].lpm_not_in_ptree_pmf = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
    }
}

__float080 * RID_Router::get_lmp(
    uint8_t iface, 
    uint8_t ptree_size) { 

    return this->lmp[iface][ptree_size].lmp;
}

__float080 RID_Router::get_lmp(
    uint8_t iface, 
    uint8_t ptree_size,
    uint8_t f) { 

    return this->lmp[iface][ptree_size].lmp[f];
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
void RID_Router::print_lmp(uint8_t iface) {

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

            _prob = this->get_lmp(iface, ptree_size, f);
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

/*
 * \brief   adds (as in 'sum', '+') a prob. value to the joint_lmp 
 *          position indexed by iface_pivots. in other words, this function does 
 *          joint_lmp[iface_pivots] += value
 *
 * multi-dimensional arrays in C/C++ are stored in row major order. what does 
 * it mean? e.g. say we have a 3D matrix, 'the_matrix' of dimension D x D x D. 
 * to access element the_matrix[1][3][5], we de-reference a pointer as such: 
 * *(the_matrix + (1 * ((D)^2)) + (3 * ((D)^1)) + (5 * ((D)^0)))
 *
 * this is why we need a function to do something as simple as 
 * joint_lmp[iface_pivots] += value
 * 
 * \param   joint_lmp               the function adds (as in 'sum', '+') value 
 *                                  to the iface_pivots position in this array.
 * \param   iface_pivots            index in joint_lmp to which the 
 *                                  prob. value should be added.
 * \param   value                   the prob. value to be added to the array.
 */
void RID_Router::add_joint_lmp_prob(
    std::map<std::string, __float080> * joint_lmp,
    int * iface_pivots,
    __float080 value) {

    // pad 1 digit numbers w/ zeros to form unique keys. e.g. [1][1][0][0][0][10]
    // the same as [1][10][0][0][1][0] without the padding
    std::ostringstream iface_pivots_str("");
    for (int i = 0; i < this->iface_num; i++) {

        if (iface_pivots[i] < 10)
            iface_pivots_str << 0;
        iface_pivots_str << iface_pivots[i];
    }

    // std::cout << "RID_Router::add_joint_lmp_prob() : [INFO] key for ";
    // for (int k = 0; k < this->iface_num; k++)
    //     std::cout << "[" << iface_pivots[k] << "]";
    // std::cout << " = " << iface_pivots_str.str() << std::endl;

    (*joint_lmp)[iface_pivots_str.str()] += value;
}

/*
 * \brief   gets a prob. value to the joint_lmp 
 *          position indexed by iface_pivots. in other words this function 
 *          returns joint_lmp[iface_pivots]
 *
 * multi-dimensional arrays in C/C++ are stored in row major order. what does 
 * it mean? e.g. say we have a 3D matrix, 'the_matrix' of dimension D x D x D. 
 * to access element the_matrix[1][3][5], we de-reference a pointer as such: 
 * *(the_matrix + (1 * ((D)^2)) + (3 * ((D)^1)) + (5 * ((D)^0)))
 *
 * this is why we need a function to do something as simple as 
 * returning joint_lmp[iface_pivots]
 * 
 * \param   joint_lmp       the function adds (as in 'sum', '+') value 
 *                                  to the iface_pivots position in this array.
 * \param   iface_pivots            index in joint_lmp to which the 
 *                                  prob. value should be added.
 * \param   value                   the prob. value to be added to the array.
 */
__float080 RID_Router::get_joint_lmp_prob(    
    std::map<std::string, __float080> * joint_lmp,
    int * iface_pivots) {

    std::ostringstream iface_pivots_str("");
    for (int i = 0; i < this->iface_num; i++) {

        if (iface_pivots[i] < 10)
            iface_pivots_str << 0;
        iface_pivots_str << iface_pivots[i];
    }

    return (*joint_lmp)[iface_pivots_str.str()];
}

/*
 * \brief compute the largest match probabilities, for each iface of the 
 * router, i.e. L_{i,p} random variables.
 *
 * L_{i,p} expresses the the prob of having a request match a forwarding entry of size 
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
int RID_Router::calc_lmp(
    uint8_t request_size, 
    int * tp_sizes, 
    __float080 * f_r_distribution) {

    // if this router already has an lmp table, don't re-calculate
    if (this->lmp_calculated) { 
        std::cout << "RID_Router::calc_lmp() : [INFO] pre-calculated lmp table available. skipping." << std::endl;
        std::cout << "RID_Router::calc_lmp() : [INFO] sizes of lmp table: "
            << "\n\t [per iface] : " << ((this->f_max + 1) * (this->f_max + 1) * sizeof(__float080)) << " byte"
            << "\n\t [TOTAL (x" << (int) this->iface_num << " ifaces)] : " 
                << (this->iface_num * (this->f_max + 1) * (this->f_max + 1) * sizeof(__float080)) << " byte" << std::endl;
        return 0;
    }

    // aux. arrays used in fp rate calculation
    __float080 * log_prob_fp_not_larger_than = (__float080 *) calloc(request_size + 1, sizeof(__float080));
    __float080 * log_prob_not_fp = (__float080 *) calloc(this->f_max, sizeof(__float080));
    // the k used in the fp rate formulas (see paper)
    __float080 k = (log(2) * ((__float080) this->bf_size)) / ((__float080) request_size);
    // variable which saves probability of not having a match of any size. 
    // used in calc_egress_ptree_probs().
    // FIXME: is this necessary? isn't the probability of a true negative 
    // accessible via other means?
    this->prob_true_negative = 1.0;

    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++) {

        // if iface isn't initialized, abort
        if (this->fwd_table[iface].iface == -1) {
            fprintf(stderr, "RID_Router::calc_lmp() : iface %d uninitialized. aborting.\n", iface);
            return -1;
        }

        // allocate memory for L_{iface,ptree_size} matrix (only for iface).
        // for coding purposes, we abbreviate L_{iface,ptree_size} as lmp.
        this->init_lmp(iface);

        // // we don't forward packets over forbidden ifaces, as such the 
        // // P(L_{forbidden iface,<for all ptree_size>} = 0) = 1.0
        // for (int ptree_size = this->f_max; ptree_size >= 0; ptree_size--) {
            
        //     std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) iface);
        //     if (it != this->no_forwarding.end())
        //         this->lmp[iface][ptree_size].lmp[0] = 1.0;
        // }

        // // skip iface if it's in the 'no forwarding' list
        // std::set<uint8_t>::iterator it = this->no_forwarding.find((uint8_t) iface);
        // if (it != this->no_forwarding.end())
        //     continue;

        // if there are no entries for this iface, a match is impossible for 
        // this [iface][ptree_size][f] combination. 
        // set [iface][ptree_size][0] = 1.0 (a 'no match' event is guaranteed) 
        // and continue to the next iface.
        if (this->fwd_table[iface].num_entries < 1) {

            for (int ptree_size = this->f_max; ptree_size >= 0; ptree_size--)
                this->lmp[iface][ptree_size].lmp[0] = 1.0;

            continue;
        }

        // pre-calculate prob of *NOT* having *ANY* FP of sizes 0 to f in iface.
        // output (log_prob_not_fp) is logarithmic.
        if (calc_log_prob_not_fp(
            (__float080) (this->bf_size), (__float080) request_size, k,
            (__float080) this->fwd_table[iface].num_entries,
            this->fwd_table[iface].f_distribution,
            f_r_distribution,
            log_prob_not_fp) < 0) {

            free(log_prob_not_fp);
            fprintf(stderr, "RID_Router::calc_lmp() : couldn't retrieve log FP rates\n");

            return -1;
        }

        // pre-calculate the probabilities of *NOT* having FPs larger than some 
        // size f. this will be useful for the L_{i,ptree_size} calculations. e.g. 
        // note that L_{iface,ptree_size} = f can only happen if a match larger 
        // than f doesn't occur.
        for (int f = (this->f_max - 1); f >= 0; f--)
            log_prob_fp_not_larger_than[f] = log_prob_fp_not_larger_than[f + 1] + log_prob_not_fp[f];

        // calculate the *LOCAL* component of the egress prefix tree 
        // probabilities
        calc_egress_ptree_probs(
            EGRESS_PTREE_PROB_MODE_LOCAL,
            iface,
            log_prob_fp_not_larger_than,
            log_prob_not_fp);

        this->prob_true_negative *= exp(log_prob_fp_not_larger_than[0]);

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
                    this->lmp[iface][ptree_size].lmp[f] = 0.0;

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
                            this->lmp[iface][ptree_size].lmp[f] = 
                                (__float080) ((this->fwd_table[iface].num_entries > 0) ? 0.0 : 1.0);

                        } else {

                            // if f >= ptree_size, the largest match will *AT 
                            // LEAST* be f.
                            this->lmp[iface][ptree_size].lmp[f] = 
                                exp(log_prob_fp_not_larger_than[f]);
                        }

                    } else {

                        // sub-case 2) no TP exists for f
                        if (f > 0) {

                            // // if f is less than the prefix tree size, then 
                            // // any match will be larger than f
                            // if (f < ptree_size) {

                            //     this->lmp[iface][ptree_size].lmp[f] = 0.0;

                            // } else if (f == ptree_size) {
                            if (f == ptree_size) {

                                // the prob of this event can be calculated as:
                                //  * not having fp matches larger than f : log_prob_fp_not_larger_than[f]
                                //  * AND having a fp match of size f : 1.0 - exp(log_prob_not_fp[f - 1])
                                this->lmp[iface][ptree_size].lmp[f] = 
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

                                    this->lmp[iface][ptree_size].lmp[f] = 0.0;
                                
                                } else {

                                    // the prob of this event can be calculated as:
                                    //  * not having fps larger than f : log_prob_fp_not_larger_than[f]
                                    //  * AND having a fp of size f : 1.0 - exp(log_prob_not_fp[f - 1])
                                    this->lmp[iface][ptree_size].lmp[f] = 
                                        exp(log_prob_fp_not_larger_than[f]) 
                                        * (1.0 - exp(log_prob_not_fp[f - 1])) 
                                        * ((this->fwd_table[iface].num_entries > 0) ? 1.0 : 0.0);
                                }
                            }

                        } else {

                            // FIXME : this also seems like something we could 
                            // cut
                            if (f < ptree_size) {

                                this->lmp[iface][ptree_size].lmp[f] = 
                                    (__float080) ((this->fwd_table[iface].num_entries > 0) ? 0.0 : 1.0);

                            } else {

                                this->lmp[iface][ptree_size].lmp[f] = 
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

    // mark the lmps as calculated
    this->lmp_calculated = true;

    return 0;
}

__float080 RID_Router::calc_log_forwarding_event_prob(
    uint8_t ptree_size,
    uint8_t ptree_iface, 
    int * iface_pivots,
    int & iterations) {

    __float080 log_joint_largest_match_prob = 0.0;
    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++) {

        if ((iface == ptree_iface))
            log_joint_largest_match_prob += log(this->get_lmp(iface, ptree_size, iface_pivots[iface]));
        else
            log_joint_largest_match_prob += log(this->get_lmp(iface, 0, iface_pivots[iface]));

        iterations++;
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
    // the individual L_{i,ptree_size} and saving it in joint_lmp.
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

            if (this->tp_sizes[i] > ptree_size)
                iface_pivots_max[i] = this->tp_sizes[i];
            else
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
    int iterations = 0;
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

                joint_prob_sum += exp(this->calc_log_forwarding_event_prob(
                    ptree_size, ptree_iface, iface_pivots, iterations));

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

            // std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] (-*) [" << curr_i << "],[" << curr_f << "]" << std::endl;

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
            // std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] (*-) [" << curr_i << "],[" << curr_f << "]" << std::endl;
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
 * \brief   compute the probability of all possible forwarding events,
 *          given the request is bound to a fp tree of size s, which continues 
 *          over iface i
 * 
 * \param   ptree_size              the joint distribution should consider this 
 *                                  'fp tree' size *ONLY*.
 * \param   ptree_iface             use the L_{i,ptree_size} values for ptree_iface, 
 *                                  L_{i,0} values for all other ifaces.
 * \param   joint_lmp               save the joint distribution computations in 
 *                                  this array.
 *
 */
int RID_Router::calc_forwarding_event_probs(
    uint8_t ptree_size,
    uint8_t ptree_iface,
    std::map<std::string, __float080> * joint_lmp) {

    // we represent joint events as a pivot array of size iface_num. e.g.
    //
    //        iface index : 0  1  2  3
    //                      v  v  v  v 
    //      iface_pivots = {0, 1, 2, 1} 
    //
    // represents the joint event in which iface 0 has no match (size 0), iface 
    // 1 has a match of size 1, iface 2 a match of size 2 and iface 3 a match 
    // of size 1. the idea is to cycle through all 
    // combinations, and calculate the probability of that event by multiplying 
    // the individual L_{i,ptree_size} and saving it in joint_lmp.
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

            if (this->tp_sizes[i] > ptree_size)
                iface_pivots_max[i] = this->tp_sizes[i];
            else
                iface_pivots_max[i] = ptree_size;

        } else {

            // ffs() gives the position of the **least** significant bit which 
            // is set. to find highest bit set, keep using ffs() and flipping 
            // the lsb set until the bitmask is 0.
            do {
                max_pivot = ffs(f_bitmask);
                f_bitmask = f_bitmask & (f_bitmask - 1);
            // } while(ffs(f_bitmask) != 0);
            } while(f_bitmask != 0);

            iface_pivots_max[i] = max_pivot;
        }
    }

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;
    int iterations[2] = {0, 0};

    __float080 prob_iface_in_ptree = 0.0;
    if (ptree_size > 0) {

        if (this->iface_on_ptree[ptree_iface]) {

            std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] ("
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

        // if the ptree_iface is def. not on ptree, then it's impossible for it 
        // to be the continuation of a prefix tree
        __float080 nr_valid_egress_ifaces = (__float080) (this->iface_num - this->no_forwarding.size());
        if (nr_valid_egress_ifaces < 1.0) nr_valid_egress_ifaces = 1.0;

        prob_iface_in_ptree = this->ingress_ptree_prob[ptree_size] * (1.0 / nr_valid_egress_ifaces);
    }

    // if prob_iface_in_ptree is 0.0, no need for further processing 
    if (prob_iface_in_ptree == 0.0)
        return 0;

    // placeholder for the individual probabilities of joint events
    __float080 log_prob = 0.0;

    // get the sum of probs of the joint events involved in this 
    // <ptree_size, ptree_iface> pair. 
    // FIXME : get_joint_prob_sum() has the same basic structure as 
    // calc_forwarding_event_probs(), and its invocation here 
    // seems overkill. unfortunately, i haven't found another way to deal 
    // with this issue.
    __float080 joint_prob_sum = calc_joint_largest_match_prob_sum(ptree_size    , ptree_iface);

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
                log_prob = this->calc_log_forwarding_event_prob(
                    ptree_size, ptree_iface, iface_pivots, iterations[1]);
                iterations[0]++;

                // scale the event probability by the following factors:
                //  1) probability of having a request coming into the router (ingress_prob)
                //  2) probability of having ptree_iface as the continuation 
                //     of a prefix tree of size ptree_size
                //  3) sum of all individual possible joint events, in order 
                //     to allow normalization of their probabilities
                log_prob += log(prob_iface_in_ptree);
                log_prob += log(this->ingress_prob);

                if (joint_prob_sum > 0.0)
                    log_prob -= log(joint_prob_sum);
                // add the calculated probability to the matrix of joint 
                // probabilities
                __float080 prob = exp(log_prob);
                this->add_joint_lmp_prob(joint_lmp, iface_pivots, prob);

                // if (prob > 0.0) {

                //     std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] JOINT_PROB" 
                //         << "(" << (int) ptree_iface << ", " << (int) ptree_size << ")";
                //     for (int k = 0; k < this->iface_num; k++)
                //         std::cout << "[" << iface_pivots[k] << "]";
                //     std::cout << " = " << prob << std::endl;

                //     std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] (SUM) JOINT_PROB";
                //     for (int k = 0; k < this->iface_num; k++)
                //         std::cout << "[" << iface_pivots[k] << "]";
                //     std::cout << " = " << this->get_joint_lmp_prob(joint_lmp, iface_pivots) << std::endl;
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

            // std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] (-*) [" << curr_i << "],[" << curr_f << "]" << std::endl;

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
            // std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] (*-) [" << curr_i << "],[" << curr_f << "]" << std::endl;
            curr_i--;
        }
    }

    std::cout << "RID_Router::calc_forwarding_event_probs() : [INFO] # of iterations : " 
        << "\n\t [INNER] : " << iterations[1]
        << "\n\t [OUTER] : " << iterations[0]
        << "\n\t [RATIO] : " << (float) iterations[1] / (float) iterations[0]
        << std::endl;

    // delete the memory allocated by the function
    free(iface_pivots);
    free(iface_pivots_max);
    free(iface_bitmasks);

    return 0;
}

/*
 * \brief   calculates marginal probability of a forwarding event of the form
 *          [*, *, f, ..., *], i.e. assuming iface i has a match of size f. this
 *          produces the probability of choosing iface i to forward a packet, 
 *          when the longest match is of size f.
 *
 *          the result also depends on a 'mode':
 *              - STRICT     : in strict longest prefix match mode,  
 *                             it computes the marginal probabilities 
 *                             over match sizes strictly less than f, for all 
 *                             ifaces other than i.
 *              - NON_STRICT : in non-strict LPM mode, it computes the 
 *                             marginal probability over match sizes less than 
 *                             or equal to f, for all ifaces. 
 *
 * 
 * \param   mode        strict or non-strict longest prefix match mode
 * \param   iface       iface which is fixed during computation of marginal prob
 * \param   iface_size  size of match on the fixed iface
 * \param   joint_lmp   hash table w/ probability of possible (i.e. w/ non-zero prob) 
 *                      forwarding event
 *
 */
__float080 RID_Router::calc_marginal_prob(
    uint8_t mode,
    uint8_t fixed_iface, 
    uint8_t fixed_iface_size,
    std::map<std::string, __float080> * joint_lmp) {

    printf("RID_Router::calc_marginal_prob() : [fixed_iface : %d][fixed_iface_size : %d][MODE : %s]\n", 
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
    } else if (this->mm_mode == MMH_RANDOM && mode == MODE_EI_INCLUSIVE) {
        distribute_probs = true;
    } else if (this->mm_mode == MMH_FALLBACK && mode == MODE_EI_EXCLUSIVE) {
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

        for (uint8_t k = 0; k < fixed_iface; k++) {

            int prev_max_pivot = 0;
            int max_pivot = 0;
            int f_bitmask = this->fwd_table[k].f_bitmask;

            do {
                prev_max_pivot = max_pivot;

                max_pivot = ffs(f_bitmask);
                if (max_pivot >= fixed_iface_size)
                    break;

                f_bitmask = f_bitmask & (f_bitmask - 1);
                
            } while(ffs(f_bitmask) != 0);

            iface_pivots_max[k] = std::min(prev_max_pivot, (int) fixed_iface_size);
        }

        for (uint8_t k = fixed_iface; k < this->iface_num; k++) {

            int fixed_iface_max_size = 0;
            int max_pivot = 0;
            int f_bitmask = this->fwd_table[k].f_bitmask;

            do {
                max_pivot = ffs(f_bitmask);
                f_bitmask = f_bitmask & (f_bitmask - 1);
            } while(ffs(f_bitmask) != 0);

            if (k == fixed_iface)
                fixed_iface_max_size = max_pivot;

            if (max_pivot > fixed_iface_max_size)
                iface_pivots_max[k] = max_pivot;
            else
                iface_pivots_max[k] = std::min(max_pivot, (int) fixed_iface_size);
        }

    } else {

        for (uint8_t k = 0; k < this->iface_num; k++) {

            int prev_max_pivot = 0;
            int max_pivot = 0;
            int f_bitmask = this->fwd_table[k].f_bitmask;

            do {
                prev_max_pivot = max_pivot;

                max_pivot = ffs(f_bitmask);
                if (max_pivot >= fixed_iface_size)
                    break;

                f_bitmask = f_bitmask & (f_bitmask - 1);

            } while(ffs(f_bitmask) != 0);

            iface_pivots_max[k] = std::min(prev_max_pivot, (int) fixed_iface_size);
        }
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
                    prob = this->get_joint_lmp_prob(joint_lmp, iface_pivots);

                    // std::cout << "RID_Router::calc_marginal_prob() : [INFO] CUMULATIVE_PROB" 
                    //     << "(" << (int) fixed_iface << ", " << (int) fixed_iface_size << ") : ";
                    // for (int k = 0; k < this->iface_num; k++)
                    //     std::cout << "[" << iface_pivots[k] << "]";
                    // std::cout << " = " << prob << std::endl;
                }

                // distribute the probabilities over the egress iface probabilities
                if (distribute_probs) {

                    int max_match = 0;
                    for (int k = 0; k < this->iface_num; k++)
                        if (iface_pivots[k] > max_match)
                            max_match = iface_pivots[k];

                    if (this->mm_mode == MMH_RANDOM) {

                        __float080 prob_multiplier = 0.0;
                        for (int k = 0; k < this->iface_num; k++) {
                            if ((iface_pivots[k]) > 0 && (iface_pivots[k] == max_match))
                                prob_multiplier += 1.0;
                        }

                        if (prob_multiplier > 0.0)
                            prob_multiplier = (1.0 / (prob_multiplier));


                        for (int k = 0; k < this->iface_num; k++) {
                            if ((iface_pivots[k]) > 0 && (iface_pivots[k] == max_match))
                                this->egress_iface_prob[k][iface_pivots[k]] += (prob_multiplier * prob);
                        }

                    } else {

                        for (int k = 0; k < this->iface_num; k++) {
                            if ((iface_pivots[k]) > 0 && (iface_pivots[k] == max_match))
                                this->egress_iface_prob[k][iface_pivots[k]] += (prob);
                        }
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

int RID_Router::calc_iface_probs(
    int ptree_size,
    std::map<std::string, __float080> * joint_lmp) {

    // // LOCAL MATCH (LLM event)
    // int local_iface_bitmask = this->fwd_table[IFACE_LOCAL].f_bitmask;
    // for (uint8_t fixed_iface_size = ffs(local_iface_bitmask); fixed_iface_size != 0; fixed_iface_size = ffs(local_iface_bitmask)) {

    //     local_iface_bitmask = local_iface_bitmask & (local_iface_bitmask - 1);
    //     this->iface_events_prob[EVENT_LLM] += this->calc_marginal_prob(
    //                                                     MODE_EI_INCLUSIVE,
    //                                                     IFACE_LOCAL,
    //                                                     fixed_iface_size,
    //                                                     joint_lmp);
    // }

    // std::cout << "RID_Router::calc_iface_probs() : [INFO]"
    //     << "\n\t P('EVENT_LLM' [ptree_size = " << (int) ptree_size << "]') = " << this->iface_events_prob[EVENT_LLM]
    //     << std::endl;

    // NO MATCHES : just call get_joint_lmp_prob() with iface_pivots = {0, 0, ..., 0}
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    this->iface_events_prob[EVENT_NLM] += this->get_joint_lmp_prob(joint_lmp, iface_pivots);
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
            prob_single = this->calc_marginal_prob(
                                    MODE_EI_EXCLUSIVE,
                                    fixed_iface,
                                    fixed_iface_size,
                                    joint_lmp);

            // in some cases, we may want the probability of having the largest 
            // match on iface, even if together with any other iface (e.g. if 
            // a flooding multiple match resolution strategy is in place). we 
            // call this MODE_EI_INCLUSIVE.
            prob_multiple = this->calc_marginal_prob(
                                    MODE_EI_INCLUSIVE,
                                    fixed_iface,
                                    fixed_iface_size,
                                    joint_lmp);

            // keep track of the total EI probability as well (sanity check)
            if (prob_single > 0.0) {

                std::cout << "RID_Router::calc_iface_probs() : [INFO]"
                    << "\n\t P('EGRESS_IFACE_EXCLUSIVE' [" << (int) fixed_iface << "][" << (int) fixed_iface_size << "]') = " << prob_single
                    << std::endl;
            }

            if (prob_multiple > 0.0) {

                std::cout << "RID_Router::calc_iface_probs() : [INFO]"
                    << "\n\t P('EGRESS_IFACE_INCLUSIVE' [" << (int) fixed_iface << "][" << (int) fixed_iface_size << "]') = " << prob_multiple
                    << std::endl;
            }

            // save the probability of a 'clean' single link match (SLM) event
            if (fixed_iface > 0)
                this->iface_events_prob[EVENT_SLM] += prob_single;
            else
                this->iface_events_prob[EVENT_LLM] += prob_multiple;

            // save the probability of multiple link match (MLM) events
            this->iface_events_prob[EVENT_MLM] += (prob_multiple - prob_single);
            std::cout << "RID_Router::calc_iface_probs() : [INFO] P('EVENT_MLM' [ptree_size = " << ptree_size << "]) = " 
                << this->iface_events_prob[EVENT_MLM] << std::endl;
            std::cout << "RID_Router::calc_iface_probs() : [INFO] P('EVENT_SLM' [ptree_size = " << ptree_size << "]) = " 
                << this->iface_events_prob[EVENT_SLM] << std::endl;
            std::cout << "RID_Router::calc_iface_probs() : [INFO] P('EVENT_LLM' [ptree_size = " << ptree_size << "]) = " 
                << this->iface_events_prob[EVENT_LLM] << std::endl;

            // // keep track of the total EI probability as well (sanity check)
            // if (prob_single > 0.0 || prob_multiple > 0.0) {

            //     std::cout << "RID_Router::calc_iface_probs() : [INFO] "
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

    return this->iface_events_prob[event];
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
void RID_Router::print_iface_events_prob() {

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

        prob = this->iface_events_prob[event];
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