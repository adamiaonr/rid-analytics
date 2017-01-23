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
    uint8_t access_tree_index, 
    uint8_t access_tree_height,
    uint8_t height, 
    uint8_t width,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max,
    uint16_t bf_size,
    int mm_mode) {

    // <tree index>.<height>.<width>
    this->id = std::to_string(access_tree_index) 
        + std::string(".") + std::to_string(height) 
        + std::string(".") + std::to_string(width);
    // other parameters
    this->access_tree_index = access_tree_index;
    this->access_tree_height = access_tree_height;
    this->height = height;
    this->width = width;
    // ???
    this->total_joint_prob = 0.0;
    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    this->f_max = f_max;
    this->bf_size = bf_size;

    this->leaf = false;

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

        for (uint8_t _f = 0; _f < this->f_max; _f++) {

            free(this->lpm_pmf[_iface][_f].lpm_pmf_prob);

            // free(this->lpm_pmf_row[_iface][_f]);
        }

        free(this->lpm_pmf[_iface]);
    }

    // free(lpm_pmf);

    free(this->iface_events_pmf);

    for (uint8_t _iface = IFACE_LOCAL + 1; _iface < this->iface_num; _iface++)
        free(this->egress_ptree_prob[_iface]);

    for (uint8_t _iface = IFACE_LOCAL + 1; _iface < this->iface_num; _iface++)
        free(this->egress_iface_prob[_iface]);
}

int RID_Router::init(
    uint8_t height, 
    uint8_t width,
    uint64_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max,
    uint16_t bf_size,
    int mm_mode) {

    if (this->initialized)
        return 0;

    // <tree index>.<height>.<width>
    this->id = std::to_string(access_tree_index) 
        + std::string(".") + std::to_string(height) 
        + std::string(".") + std::to_string(width);
    // other parameters
    this->height = height;
    this->width = width;
    // ???
    this->total_joint_prob = 0.0;
    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    this->f_max = f_max;
    this->bf_size = bf_size;

    this->leaf = false;

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

int RID_Router::forward(
    uint8_t request_size,
    uint8_t ingress_iface,
    int * tp_sizes,
    __float080 ingress_prob,
    __float080 * ingress_ptree_prob,
    __float080 * f_r_distribution) {

    // clear the no forwarding list 
    this->no_forwarding.clear();
    // take note of the ingress iface for further use. bottom line: we don't 
    // want to forward over it again and create a loop.
    this->no_forwarding.insert(ingress_iface);

    // special case happens with peering relationships: you don't want to forward 
    // over an upstream iface and risk creating a loop.
    std::cout << "RID_Router::forward() : [INFO] forwarding exclusions: " 
        << "\n\t" << this->get_id() << "[" << (int) ingress_iface << "]@[" << (int) this->get_height() << "] -> " 
        << this->get_next_hop(ingress_iface).router->get_id() 
            << "@[" << (int) this->get_next_hop(ingress_iface).router->get_height() << "]" << std::endl;

    if ((this->get_next_hop(ingress_iface).router->get_height() <= this->height)
        && (this->get_id() != "0.3.0")) {

        // we never consider IFACE_LOCAL for exclusion though...
        for (int i = (IFACE_LOCAL + 1); i < this->iface_num; i++) {
            // FIXME: this will add the ingress iface again
            if (this->get_next_hop(i).router->get_height() <= this->height)
                this->no_forwarding.insert(i);
        }
    }

    // also, if the iface doesn't have any entries, add it to the 'no flight' 
    // list
    for (int i = IFACE_LOCAL; i < this->iface_num; i++)
        if (this->fwd_table[i].num_entries == 0)
            this->no_forwarding.insert(i);

    std::cout << "RID_Router::forward() : [INFO] not forwarding over ifaces {";
    for (std::set<int>::iterator it = this->no_forwarding.begin();
        it != this->no_forwarding.end();
        ++it) {

        std::cout << (int) (*it) << ", ";        
    }
    std::cout << "}" << std::endl;

    // save the ingress probabilities for each diff. prefix tree size (i.e. 
    // the probability that a packet is already bound to a tree of size |F|)
    this->ingress_prob = ingress_prob;
    std::cout << "RID_Router::forward() : INGRESS_PTREE_PROB = " 
        << this->ingress_prob << std::endl;

    for (uint8_t ptree_size = 0; ptree_size <= this->f_max; ptree_size++) {
        this->ingress_ptree_prob[ptree_size] = ingress_ptree_prob[ptree_size];

        std::cout << "RID_Router::forward() : INGRESS_PTREE_PROB[" 
            << (int) ptree_size << "] = " << this->ingress_ptree_prob[ptree_size] << std::endl;
    }
    
    // 1) calculate the distribution for the L_{i,p} random variables, for 
    // each iface i and 'prefix tree' size p. 
    // L_{i,p} expresses the the prob of having 
    // a match of size L, at iface i, considering a binding to a prefix 
    // tree of size p.
    if (calc_lpm_pmf(
        request_size, 
        tp_sizes,              // P('prefix tree size = p')
        f_r_distribution) < 0) {

        return -1;
    }

    // print L_{i,p} for all (i,p) pairs (ingress_iface should be 0)
    for (uint8_t iface = IFACE_LOCAL; iface < this->iface_num; iface++)
        this->print_lpm_pmf(iface);

    std::cout << "RID_Router::forward() : egress ptree size probs:" << std::endl;
    for (int i = 0; i < this->iface_num; i++) {
        for (int f = 0; f <= this->f_max; f++) {

            std::cout << "\tEGRESS_PTREE_PROB[" << (int) i << "][" << (int) f << "] = " 
                << this->egress_ptree_prob[i][f] << std::endl;
        }
    }

    // 2) calculate the joint probability distribution of the 
    // random variables L_{0,p} x L_{1,p} x ... x L_{|R|_{max},p}.

    // 2.1) initialize a single array to compute the joint prob distribution for 
    // all 'prefix tree' sizes ptree_size
    __float080 * joint_lpm_matrix = NULL;
    init_joint_lpm_pmf(&(joint_lpm_matrix));

    // initialize iface_events_pmf
    for (int i = 0; i < EVENT_NUM; i++)
        this->iface_events_pmf[i] = 0.0;
    // initialize egress iface probs
    for (int i = 0; i < this->iface_num; i++)
        for (int f = 0; f < (this->f_max + 1); f++)
            this->egress_iface_prob[i][f] = 0.0;
    // reset the tracker of added joint events
    this->added_pivots.clear();

    // 2.3) compute the joint distr. for p, considering each diff. iface 
    // as the iface associated with the prefix tree of a certain size
    for (int ptree_size = 0; ptree_size <= this->f_max; ptree_size++) {

        // if the probability of having a request bound to a ptree of size 
        // ptree_size is 0.0, then it is impossible for an upstream 
        // iface to be 'bound' to this ptree_size. thus, we skip this 
        // calculation.
        if (ingress_ptree_prob[ptree_size] == 0.0) {
            std::cout << "RID_Router::forward() : skipping ptree[" 
                << ptree_size << "]" << std::endl;
            continue;
        }

        // re-use the joint probability matrix for the computation of 
        // joint prob of each prefix tree size
        clear_joint_lpm_pmf(&(joint_lpm_matrix));

        for (uint8_t ptree_iface = 0; ptree_iface < this->iface_num; ptree_iface++) {

            // we assume there are no loops, and as such the 'fake' prefix 
            // tree will not be associated with the ingress iface. notice that 
            // if the packet is indeed bound to a 'fake' prefix tree, the 
            // corresponding FP entry should be announced by an upstream router 
            // somewhere
            std::set<int>::iterator it = this->no_forwarding.find((int) ptree_iface);
            if (it != this->no_forwarding.end())
                continue;

            // we'll never have a non-leaf router to be the origin of a 
            // prefix tree, so skip the local iface.
            // FIXME: this might be a source of problems...
            if ((ptree_iface == IFACE_LOCAL) 
                && ((this->height < this->access_tree_height)
                    || (this->id == "0.3.0")))
                continue;

            calc_joint_lpm_pmf(joint_lpm_matrix, ptree_size, ptree_iface);
        }

        // 3) calc the prob distribution for ifaces, for each ptree_size
        if (calc_iface_events_pmf(joint_lpm_matrix, ptree_size) < 0)
            return -1;
    }

    // we now deduct the probability of exclusive iface matches to correctly 
    // calculate the probability of multiple matches
//    this->iface_events_pmf[EVENT_MLM] -= this->iface_events_pmf[EVENT_SLM];

    std::cout << "RID_Router::forward() :"
        << "\n\t P('SIM') = " << this->iface_events_pmf[EVENT_SLM]
        << "\n\t P('MIM') = " << this->iface_events_pmf[EVENT_MLM]
        << std::endl;

    // FIXME : i don't know if this is the correct way of doing this...
    // if the request is going up, we relay NIS cases to the router upstream. 
    // once we reach the top of the topology, we give up relaying and use 
    // a 'fallback' to origin
    if ((this->height > 0) &&
        ((this->get_id() == "0.3.0") || (this->get_next_hop(ingress_iface).router->get_height() > this->height))) {

        this->egress_iface_prob[1][0] += this->iface_events_pmf[EVENT_NLM];

        __float080 egress_ptree_prob_sum = 0.0;
        if ((this->get_id() != "0.3.0")) {

            this->egress_ptree_prob[1][0] += this->iface_events_pmf[EVENT_NLM];

                if (this->egress_ptree_prob[1][0] > 0.0) {

                    for (int f = 0; f <= this->f_max; f++)
                        egress_ptree_prob_sum += this->egress_ptree_prob[1][f];

                    for (int f = 0; f <= this->f_max; f++)
                        this->egress_ptree_prob[1][f] = (this->egress_ptree_prob[1][f] / egress_ptree_prob_sum);
            }

            int e = this->no_forwarding.erase(1);
            std::cout << "RID_Router::forward() : [INFO] erased " << e << " elements from blocked ifaces" << std::endl;
        }
    }

    std::cout << "RID_Router::forward() : [INFO] P('EVENT_MLM' / valid_ifaces) = " 
        << this->iface_events_pmf[EVENT_MLM] / 3.0 << std::endl;

    this->print_iface_events_pmf();
    this->print_egress_iface_prob();

    // think it's safe to do this? think again...
    free(joint_lpm_matrix);

    return 0;
}

int RID_Router::add_fwd_table_entry(
    int iface, 
    __float080 iface_proportion, 
    std::map<int, __float080> size_dist,
    RID_Router * next_hop_router,
    int next_hop_iface) {

    if (iface < (int) IFACE_LOCAL || iface > (this->iface_num - 1)) {
        std::cerr << "RID_Router::add_fwd_table_entry() : [ERROR] invalid iface index: "
            << iface << std::endl;
        return -1;
    }

    // the nr. of entries for this iface is a fraction of the fwd table size
    this->fwd_table[iface].num_entries = (__float080) this->fwd_table_size * iface_proportion;    
    this->fwd_table[iface].iface = iface;
    this->fwd_table[iface].iface_proportion = iface_proportion;    

    this->fwd_table[iface].f_distribution = 
        (__float080 *) calloc(this->f_max, sizeof(__float080));

    for (int f = 0; f < this->f_max; f++)
        this->fwd_table[iface].f_distribution[f] = size_dist[f + 1];

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

int RID_Router::get_num_entries(uint8_t iface) {
    return (int) this->fwd_table[iface].num_entries;
}

uint64_t RID_Router::get_fwd_table_size() { 
    return this->fwd_table_size; 
}

uint8_t RID_Router::get_iface_num() { 
    return this->iface_num; 
}

uint8_t RID_Router::get_f_max() { 
    return this->f_max; 
}

uint8_t RID_Router::get_access_tree_index() { 
    return this->access_tree_index; 
}

uint8_t RID_Router::get_height() { 
    return this->height; 
}

uint8_t RID_Router::get_width() { 
    return this->width; 
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

int RID_Router::get_log_fp_rates(
    __float080 m, 
    __float080 n, 
    __float080 k,
    __float080 f_entries,               // nr. of entries for iface
    __float080 * f_distribution, 
    __float080 * f_r_distribution,
    __float080 * log_fp_rates) {        // function fills this array

    __float080 _n_entries = 0.0;
    __float080 _f_entries = 0.0; 
    __float080 _subtotal_f_r = 0.0; 

    // as a precaution, reset log_fp_rates to 0.0
    for (int f = 0; f < this->f_max; f++)
        log_fp_rates[f] = 0.0;

    for (int _f = 0; _f < this->f_max; _f++) {

        _f_entries = f_entries * f_distribution[_f];
        // std::cout << "RID_Router::get_log_fp_rates() : _f_entries[" << _f << "] = " << _f_entries << " (out of " << f_entries << ")" << std::endl;

        if (_f == 0) {

            log_fp_rates[_f] += _f_entries * log(1.0 - fp_rate(m, n, k, (__float080) (_f + 1)));  

        } else {

            _subtotal_f_r = 0.0;

            for (int _jf = 0; _jf <= _f; _jf++) {
                _subtotal_f_r += f_r_distribution[_jf];
            }

            for (int _jf = 0; _jf <= _f; _jf++) {

                _n_entries = (__float080) _f_entries * (f_r_distribution[_jf] / _subtotal_f_r);
                log_fp_rates[_f] += _n_entries * log(1.0 - fp_rate(m, n, k, (__float080) (_jf + 1)));

                // std::cout << "RID_Router::get_log_fp_rates() : log_fp_rates[" 
                //     << _f << "][" << _jf << "] = " 
                //     << exp(_n_entries * log(1.0 - fp_rate(m, n, k, (__float080) (_jf + 1)))) << std::endl;

                // std::cout << "RID_Router::get_log_fp_rates() : log_fp_rates[" 
                //     << _f << "][" << _jf << "] = " 
                //     << exp(10 * log(1.0 - fp_rate(m, n, k, (__float080) (_jf + 1)))) << std::endl;
            }
        }

        // std::cout << "RID_Router::get_log_fp_rates() : log_fp_rates[" 
        //     << _f << "] = " << exp(log_fp_rates[_f]) << std::endl;
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

/*
 * \brief   allocates memory for an array capable of holding a joint distribution 
 *          matrix for all L_i,p RVs (i.e. for all ifaces in the router)
 *
 * since n-dimensional arrays in C/C++ are initialized in row major order, we 
 * access the contents of joint_f_pmf by 'smartly' dereferencing the pointer 
 * (wooo, watch out, we've got a rocket scientist over here...)
 * 
 * \param   joint_prob_matrix   the array to initialize
 */
void RID_Router::init_joint_lpm_pmf(
    __float080 ** joint_prob_matrix) {

    *joint_prob_matrix = (__float080 *) calloc(pow((this->f_max + 1), this->iface_num), sizeof(__float080));
}

/*
 * \brief   sets the joint prob distribution array to 0.0s. this is just a 
 *          wrapper for memset().
 * 
 * \param   joint_prob_matrix   the array to clear
 */
void RID_Router::clear_joint_lpm_pmf(__float080 ** joint_prob_matrix) {

    this->total_joint_prob = 0.0;
    memset(*joint_prob_matrix, 0, pow((this->f_max + 1), this->iface_num) * sizeof(__float080));
}

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
    __float080 * joint_prob_matrix,
    int * iface_pivots,
    __float080 value) {

    uint32_t joint_prob_add = 0;

    for (uint8_t i = 0; i < this->iface_num; i++) {

        joint_prob_add += iface_pivots[i] * pow((this->f_max + 1), i);
    }

    *(joint_prob_matrix + joint_prob_add) += value;
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
    __float080 * joint_prob_matrix,
    int * iface_pivots) {

    uint32_t joint_prob_add = 0;

    for (int i = 0; i < this->iface_num; i++) {

        joint_prob_add += iface_pivots[i] * pow((this->f_max + 1), i);
    }

    return *(joint_prob_matrix + joint_prob_add);
}

/*
 * \brief   compute the distribution of L_{i,p} random variables, for all 
 *          <iface, prefix tree size> pairs in the router. L_{i,p} doesn't 
 *          take true positive info into account.
 *
 * the random variable L_{i,p} represents the longest match reported by an 
 * iface i, assuming a prefix tree of size p. e.g., P(L_{0,1} = 5) is the prob 
 * of getting a longest match equal to 5, for iface 0 and assuming that the 
 * request got 'stuck' to a prefix tree of size 1 in some previous point 
 * in the path.
 *
 * the objective is to build a table like the following, one for each iface:
 *
 *              |   size of longest match for iface i    |
 *              ------------------------------------------
 * | ptree size | |F| = 0 (no match) | 1 | 2 | 3 | 4 | 5 |
 * -------------------------------------------------------
 * |          0 |         "          | " | " | " | " | " |
 * |          1 |         "          | " | " | " | " | " |
 * |          2 |         "          | " | " | " | " | " |
 * |          3 |         "          | " | " | " | " | " |
 * |          4 |         "          | " | " | " | " | " |
 * |          5 |         "          | " | " | " | " | " |
 *
 *
 * L_{i,p} as if true positives don't exist. we apply true positive info later.
 *  
 * for example, P(L_{0,1} = 3) can be > 0.0, when iface 0 has a true 
 * positive entry of size 4. in rigor, P(L_{0,1} = 3) should be 0.0 because 
 * we know for sure that the true positive entry will always be a match, and 
 * thus the size of the matches for iface 0 will ALWAYS be larger than or equal 
 * to 4. HOWEVER, we first calculate L_{i,p} as if true positives don't exist, 
 * and then apply true positive information later.
 *
 * why do we do things this way?  
 * 
 * \param   fp_size_prob    probabilities of the FP match sizes that led the 
 *                          request to this router.
 *
 * \return  1 if bit i is set, 0 if not set.
 */
int RID_Router::calc_lpm_pmf(
    uint8_t request_size, 
    int * tp_sizes, 
    __float080 * f_r_distribution) {

    // probability of having an iface as part of a 'prefix tree' of size |F| = f

    // speed things up by keeping the (~FP5 * ~FP4 * ~FP4 * ... * ~FPi)
    // probabilities in an array
    __float080 * _log_prob_not_fp = (__float080 *) calloc(request_size + 1, sizeof(__float080));
    __float080 * _log_fp_rates = (__float080 *) calloc(this->f_max, sizeof(__float080));

    __float080 _k = (log(2) * ((__float080) this->bf_size)) / ((__float080) request_size);

    __float080 egress_ptree_prob_sum = 0.0;

    // FIXME: the iface is associated w/ a prefix tree at random, w/ 
    // probability 1 / <nr. of valid egress ifaces>
    int valid_ifaces = (int) this->iface_num - this->no_forwarding.size();
    if (valid_ifaces < 1)
        valid_ifaces = 1;

    // starting with _iface = IFACE_LOCAL (i.e. cache)
    for (uint8_t _iface = IFACE_LOCAL; _iface < this->iface_num; _iface++) {

        // if _iface isn't initialized, abort
        if (this->fwd_table[_iface].iface == -1) {

            fprintf(stderr, "RID_Router::calc_lpm_pmf() : iface %d uninitialized. aborting.\n", _iface);
            return -1;
        }

        // allocate memory for L_{_iface,p} matrix (only for _iface)
        this->init_lpm_pmf(_iface);

        for (int _ptree_size = this->f_max; _ptree_size >= 0; _ptree_size--) {
            
            // do not forward packets over forbidden ifaces
            std::set<int>::iterator it = this->no_forwarding.find((int) _iface);
            if (it != this->no_forwarding.end())
                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[0] = 1.0;

            // // FIXME: as a special fixme, if the next hop is 'NULL', treat it 
            // // in the same way as the ingress iface
            // if (_iface > IFACE_LOCAL && this->get_next_hop(_iface).router == NULL) {

            //     this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[0] = 1.0;
            // }
        }

        // thus, skip everything else if _iface is in the 'no flight' list
        // do not forward packets over forbidden ifaces
        std::set<int>::iterator it = this->no_forwarding.find((int) _iface);
        if (it != this->no_forwarding.end())
            continue;
        // if (_iface > IFACE_LOCAL && this->get_next_hop(_iface).router == NULL) 
        //     continue;

        // get the 1 x |R|_max matrix with the probs of *NOT* having FPs among  
        // entries of size f in _iface. these probs are given in log() 
        // form, since we'll be dealing with small numbers.

        if (get_log_fp_rates(
            (__float080) (this->bf_size), (__float080) request_size, _k,   // BF parameters
            (__float080) this->fwd_table[_iface].num_entries,   // # of entries for _iface
            this->fwd_table[_iface].f_distribution,     // distr. of diff. |F| for _iface
            f_r_distribution,                           // distr. of |F\R| for _iface
            _log_fp_rates) < 0) {                       // array were result is stored

            // FAIL! so free memory of log_fp_rates 
            free(_log_fp_rates);
            fprintf(stderr, "RID_Router::calc_lpm_pmf() : couldn't retrieve log FP rates\n");

            return -1;
        }

        // fill the (~FP5 * ~FP4 * ~FP4 * ... * ~FPi) array (the ancient 
        // technique of 'memoization', woooo...)

        // these are the probabilities of *NOT* having FPs for |F| > _f. this 
        // is present in all sub-cases (TPs, local FPs, getting stuck in 
        // 'prefix trees'). this makes sense: L_{_iface,p} = _f can only happen 
        // if a match larger than _f doesn't happen.

        // for (int f = (this->f_max - 1); f >= 0; f--) {

        //     std::cout << "RID_Router::calc_lpm_pmf() : [INFO] P(|FP| == f)[" 
        //         << (int) _iface << "][" << f << "] = " << _log_fp_rates[f]
        //         << std::endl;            
        // }

        for (int _f = (this->f_max - 1); _f >= 0; _f--) {
            _log_prob_not_fp[_f] = _log_prob_not_fp[_f + 1] + _log_fp_rates[_f];

            // std::cout << "RID_Router::calc_lpm_pmf() : [INFO] P(~|FP| > f)[" 
            //     << (int) _iface << "][" << _f << "] = " << _log_prob_not_fp[_f]
            //     << std::endl;
        }

        // we now calculate the egress prefix tree probabilities, i.e. the 
        // probability of having a packet leave a router over iface i AND in 
        // a 'wrong' prefix tree of size p. 
        // this only accounts w/ the effect of 
        // FP matches, and does not take TPs into account. this can happen 
        // if 1 of 2 events (a or b) is verified:
        //
        //  a.1) the prefix tree of size p is associated w/ iface i 
        //  a.2) a local FP match of size f \in {0, 1, ..., p} occurs
        //
        //  b.1) iface i is associated w/ a prefix tree of size q < p
        //  b.2) a local FP match of size f == p occurs
        //
        // FIXME: this could be a source of trouble...

        // FIXME: initialize the egress FP tree probs
        for (int f = 0; f <= this->f_max; f++)
            this->egress_ptree_prob[_iface][f] = 0.0;

        for (int f = 0; f <= this->f_max; f++) {

            __float080 aux = 1.0;
            if (f > 0)
                aux = ((this->fwd_table[_iface].f_distribution[f - 1] > 0.0) ? 1.0 : 0.0);

            this->egress_ptree_prob[_iface][f] += 
                this->ingress_ptree_prob[f]     // prob of _iface in prefix tree of size f
                * exp(_log_prob_not_fp[0])      // not having a FP > 0 (i.e. not having FPs at all)
                * ((this->fwd_table[_iface].num_entries > 0) ? 1.0 : 0.0)   // making sure that there are entries in this iface
                * aux;

            // if (this->ingress_ptree_prob[f] > 0.0) {

            //     if (exp(_log_prob_not_fp[0]) > 0.0)
            //         std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (1) EP["
            //             << (int) _iface << "][" << f << "] _log_prob_not_fp[0] = " << exp(_log_prob_not_fp[0]) << std::endl;

            //     if (exp(_log_prob_not_fp[0]) > 0.0)
            //         std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (1) EP["
            //             << (int) _iface << "][" << f << "] f_distribution[" << (int) _iface << "][" << f << "] = " << this->fwd_table[_iface].f_distribution[f - 1] << std::endl;

            //     if (exp(_log_prob_not_fp[0]) > 0.0)
            //         std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (1) EP["
            //             << (int) _iface << "][" << f << "] = " << this->egress_ptree_prob[_iface][f] << std::endl;
            // }
        }

        for (int f = 1; f <= this->f_max; f++) {
            for (int p = f; p >= 0; p--) {

                if (p < f) {

                    this->egress_ptree_prob[_iface][f] +=
                        this->ingress_ptree_prob[p]     // prob of _iface in prefix tree of size p
                        * exp(_log_prob_not_fp[f])      // not having a FP > f
                        * (1.0 - exp(_log_fp_rates[f - 1]))                         // generating a FP of size f
                        * ((this->fwd_table[_iface].num_entries > 0) ? 1.0 : 0.0)   // making sure that there are entries in this iface
                        * ((this->fwd_table[_iface].f_distribution[f - 1] > 0.0) ? 1.0 : 0.0);

                    // if (this->ingress_ptree_prob[p] > 0.0) {

                    //     if (exp(_log_prob_not_fp[f]) > 0.0) {

                    //         if((1.0 - exp(_log_fp_rates[f - 1])) > 0.0) {             

                    //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (2) EP["
                    //                 << (int) _iface << "][" << f << "] _log_prob_not_fp[" << f << "] = " << exp(_log_prob_not_fp[f]) << std::endl;

                    //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (2) EP["
                    //                 << (int) _iface << "][" << f << "] 1.0 - exp(_log_fp_rates[" << f - 1 << "]) = " << (1.0 - exp(_log_fp_rates[f - 1])) << std::endl;

                    //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (2) EP["
                    //                 << (int) _iface << "][" << f << "] f_distribution[" << (int) _iface << "][" << f << "] = " << this->fwd_table[_iface].f_distribution[f - 1] << std::endl;

                    //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (2) EP["
                    //                 << (int) _iface << "][" << f << "] = " << this->egress_ptree_prob[_iface][f] << std::endl;
                    //         }
                    //     }
                    // }

                } else {

                    for (int g = p; g >= 1; g--) {

                        this->egress_ptree_prob[_iface][f] +=
                            this->ingress_ptree_prob[p]             // prob of _iface in prefix tree of size p
                            * exp(_log_prob_not_fp[g])              // not having a FP > f
                            * (1.0 - exp(_log_fp_rates[g - 1]))     // generating a FP of size f
                            * ((this->fwd_table[_iface].num_entries > 0) ? 1.0 : 0.0)   // making sure that there are entries in this iface
                            * ((this->fwd_table[_iface].f_distribution[f - 1] > 0.0) ? 1.0 : 0.0);

                        // if (this->ingress_ptree_prob[p] > 0.0) {

                        //     if (exp(_log_prob_not_fp[g]) > 0.0) {

                        //         if((1.0 - exp(_log_fp_rates[g - 1])) > 0.0) {      

                        //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (3) EP["
                        //                 << (int) _iface << "][" << f << "] _log_prob_not_fp[" << g << "] = " << exp(_log_prob_not_fp[g]) << std::endl;

                        //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (3) EP["
                        //                 << (int) _iface << "][" << f << "] 1.0 - exp(_log_fp_rates[" << g - 1 << "]) = " << (1.0 - exp(_log_fp_rates[g - 1])) << std::endl;

                        //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (3) EP["
                        //                 << (int) _iface << "][" << f << "] f_distribution[" << (int) _iface << "][" << g << "] = " << this->fwd_table[_iface].f_distribution[g - 1] << std::endl;

                        //             std::cout << "RID_Router::calc_lpm_pmf() : [INFO] (3) EP["
                        //                 << (int) _iface << "][" << f << "] = " << this->egress_ptree_prob[_iface][f] << std::endl;
                        //         }
                        //     }
                        // }
                    }
                }
            }
        }

        for (int f = 0; f <= this->f_max; f++)
            egress_ptree_prob_sum += this->egress_ptree_prob[_iface][f];

        if (egress_ptree_prob_sum > 0.0)
            for (int f = 0; f <= this->f_max; f++)
                this->egress_ptree_prob[_iface][f] = (this->egress_ptree_prob[_iface][f] / egress_ptree_prob_sum);

        egress_ptree_prob_sum = 0.0;


        // if (this->fwd_table[_iface].num_entries > 0) {
        //     this->egress_ptree_prob[_iface][0] *= exp(_log_prob_not_fp[0]);
        // } else {
        //     this->egress_ptree_prob[_iface][0] = 1.0;
        // }

        // start with P(L_{_iface,_ptree_size} = |R|_max AND 
        // _ptree_size = |R|_max), go down from there...
        for (int _ptree_size = this->f_max; _ptree_size >= 0; _ptree_size--) {

            for (int _f = this->f_max; _f >= 0; _f--) {

                // FIXME : if there are no entries for this iface (it could 
                // happen given appropriate distributions), a match is impossible 
                // for this [_iface][_ptree_size][_f] combination. don't forget 
                // to set [_iface][_ptree_size][0] = 1.0 (a 'no match' event 
                // is guaranteed)
                if (this->fwd_table[_iface].num_entries == 0) {

                    if (_f > 0) {

                        this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 0.0;

                    } else {

                        this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 1.0;
                    }

                    continue;
                }

                // if _f is smaller than the max. true positive match 
                // (tp_sizes[_iface]) then _f *WON'T BE* the longest match for 
                // sure : notice we're *GUARANTEED* to have a match of at least 
                // size tp_sizes[_iface]. 
                if (_f < tp_sizes[_iface]) {

                    // _f should also be smaller than the current _ptree_size, 
                    // since _f >= _ptree_size would guarantee a positive 
                    // match of *AT LEAST* _ptree_size
                    if (_f < _ptree_size)
                        this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 0.0;

                    // FIXME: you may need to do something about the egress 
                    // probabilities here...

                } else {

                    // division of sub-cases (deterministic):
                    //  1) TP exists for |F| = _f
                    //  2) no TP exists for |F| = _f
                    if (_f == tp_sizes[_iface]) {

                        // 1) TP exists for |F| = _f

                        // if _iface is in the ptree and _f < _ptree_size : 
                        // since _ptree_size < tp_sizes[_iface], _f will never 
                        // be the longest match
                        if (_f < _ptree_size) {

                            // FIXME: this seems hack-ish : here i check if 
                            // _iface is associated with any entries
                            this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                (__float080) ((this->fwd_table[_iface].num_entries > 0) ? 0.0 : 1.0);

                        } else {

                            // if _f >= _ptree_size, the longest match will *AT 
                            // LEAST* be _f.
                            this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                exp(_log_prob_not_fp[_f]);
                        }

                    } else {

                        // 2) no TP exists for |F| = _f

                        if (_f > 0) {

                            if (_f < _ptree_size) {

                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 0.0;

                            } else if (_f == _ptree_size) {

                                // in this case, a longest match of size f 
                                // happens if:
                                //  1) no FPs w/ size > f are spontaneously 
                                //     generated
                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    (__float080) exp(_log_prob_not_fp[_f]) * ((this->fwd_table[_iface].num_entries > 0) ? 1.0 : 0.0);

                                // //  2) if the iface is in the prefix tree of 
                                // //     size _ptree_size
                                // this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] =
                                //     this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f]
                                //     * (1.0 / ((__float080) valid_ifaces)) * this->ingress_ptree_prob[_ptree_size];

                            } else {

                                // if f is larger than the current prefix tree 
                                // size, we could have a longest match of size 
                                // f iif a false positive is spontaneously 
                                // generated.

                                // the prob of this event can be calculated as:
                                //  * not having fps larger than f : _log_prob_not_fp[_f]
                                //  * AND having a fp of size f : 1.0 - - exp(_log_fp_rates[_f - 1])
                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    exp(_log_prob_not_fp[_f]) * (1.0 - exp(_log_fp_rates[_f - 1])) * ((this->fwd_table[_iface].num_entries > 0) ? 1.0 : 0.0);

                                // this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] =
                                //     this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f]
                                //     * (1.0 / ((__float080) valid_ifaces)) * this->ingress_ptree_prob[_ptree_size];
                            }

                        } else {

                            if (_f < _ptree_size) {

                                // FIXME: this seems hack-ish : here i check if 
                                // _iface is associated with any entries. if it 
                                // is, then
                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    (__float080) ((this->fwd_table[_iface].num_entries > 0) ? 0.0 : 1.0);

                            } else {

                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    exp(_log_prob_not_fp[_f]);

                                // this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] =
                                //     this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f]
                                //     * (1.0 / ((__float080) valid_ifaces)) * this->ingress_ptree_prob[_ptree_size];
                            }
                        }
                    }
                }
            }
        }
    }

    free(_log_prob_not_fp);
    free(_log_fp_rates);

    return 0;
}

void RID_Router::set_leaf() {
    this->leaf = true;
}

bool RID_Router::is_leaf() {
    return this->leaf;
}

/*
 * \brief   compute the joint probability distribution of L_{i,p} x 
 *          L_{j,p} x ... x L_{n,p} RVs, between all ifaces in the router and 
 *          for a 'prefix tree' size p.
 * 
 * \param   joint_prob_matrix       save the joint distribution computations in 
 *                                  this array. this is a value-result arg in 
 *                                  a way.
 * \param   ptree_size              the joint distribution should consider this 
 *                                  'prefix tree' size *ONLY*.
 * \param   ptree_iface             use the L_{i,ptree_size} values for ptree_iface, 
 *                                  L_{i,0} values for all other ifaces.
 *
 * \return  0 if calculation is successful, an error code < 0 if errors occur.
 */
int RID_Router::calc_joint_lpm_pmf(
    __float080 * joint_prob_matrix,
    uint8_t ptree_size,
    uint8_t ptree_iface) {

    // hold the pivots of the different iface arrays 
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int)); 

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;

    __float080 _log_prob = 0.0;

    // we use the f_distributions to know at which ifaces ptrees of size 
    // ptree_size can continue. e.g. if iface.f_distribution[ptree_size] = 0.0, 
    // meaning that iface doesn't have any entries of size ptree_size 
    // assigned to it, then it's not possible for a ptree to continue through 
    // it.
    __float080 f_dist = 0.0;
    if (ptree_size > 0)
        f_dist = this->fwd_table[ptree_iface].f_distribution[ptree_size - 1];
    else
        f_dist = 1.0;

    while (curr_i > -1) {

        // the inner loop is where all the work gets done
        if (curr_i == (this->iface_num - 1)) {

            for (int i = 0; i < (this->f_max + 1); i++) {

                // update the |F| size pivot for iface i 
                iface_pivots[this->iface_num - 1] = i;

                // calculate the log probability for the joint event 
                // encoded in iface_pivots
                _log_prob = this->calc_joint_log_prob(ptree_size, ptree_iface, iface_pivots);
                // _log_prob is scaled by P("having the 'prefix tree' of size 
                // _ptree_size assigned to _iface")

                // 2 ways for a ptree of size ptree_size to continue on iface:
                //  1) iface.f_dist[ptree_size] > 0.0
                //  2) OR the router is a leaf node (if the packet traveled 
                //     all the way to a leaf node, then its delivery is mandatory)
                if (f_dist > 0.0 || this->is_leaf())
                    _log_prob = exp(_log_prob) * this->ingress_ptree_prob[ptree_size];
                else
                    _log_prob = 0.0;

                // set this value in the joint_lpm_size_matrix 
                // note that this function adds exp(log_prob) to the value 
                // already held by joint_lpm_size_matrix in the position indexed 
                // by iface_pivots
                this->add_joint_lpm_prob(joint_prob_matrix, iface_pivots, _log_prob);
                this->total_joint_prob += _log_prob;

                if (_log_prob > 0.0) {

                    std::cout << "RID_Router::calc_joint_lpm_pmf() : [INFO] JOINT_PROB" 
                        << "(" << (int) ptree_iface << ", " << (int) ptree_size << ")";
                    for (int k = 0; k < this->iface_num; k++)
                        std::cout << "[" << iface_pivots[k] << "]";
                    std::cout << " = " << _log_prob << std::endl;

                    std::cout << "RID_Router::calc_joint_lpm_pmf() : [INFO] (SUM) JOINT_PROB";
                    for (int k = 0; k < this->iface_num; k++)
                        std::cout << "[" << iface_pivots[k] << "]";
                    std::cout << " = " << this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots) << std::endl;
                }

                std::set<int>::iterator it = this->no_forwarding.find((int) curr_i);
                if (it != this->no_forwarding.end()) {
                    iface_pivots[curr_i] = this->f_max;
                    break;
                }
            }
        }

        // if an iface is in the 'no forwarding' list, we don't need to 
        // calculate probabilities pivot values other than 0: we already 
        // know these are gonna be equal to 0.0.
        // as such, if curr_i is in no_forwarding AND we're NOT re-setting 
        // it's pivot value (i.e. setting iface_pivots[curr_i] = 0 after 
        // ++curr_f), we force its pivot to f_max, so that we skip all 
        // possible pivot values between 1 and f_max. 
        std::set<int>::iterator it = this->no_forwarding.find((int) curr_i);
        if ((it != this->no_forwarding.end()) && (iface_pivots[curr_i] > -1))
            iface_pivots[curr_i] = this->f_max;

        // extract the current pivot value for the curr_i level
        curr_f = iface_pivots[curr_i];

        if (++curr_f < (this->f_max + 1)) {

            // we're at some level curr_i, and we need to increment 
            // it's pivot counter (to curr_f)
            iface_pivots[curr_i] = curr_f;

            // after updating iface_pivots[curr_i], go back down 
            // to the level below.
            curr_i++;

            // std::cout << "RID_Router::calc_joint_lpm_pmf() : [INFO] [CURR_I][CURR_F] = " 
            //     << "[" << curr_i << "],[" << curr_f << "]" << std::endl;

        } else {

            // we've completed the iterations for the curr_i level. we 
            // now: 
            //  * reset iface_pivots[curr_i] to 0
            curr_f = -1;
            iface_pivots[curr_i] = curr_f;
            //  * go up one level (i.e. to curr_i - 1) to increment 
            //      iface_pivots[curr_i - 1]
            curr_i--;
        }
    }

    // delete the memory allocated by the function
    free(iface_pivots);

    return 0;
}

int RID_Router::get_valid_ifaces(int ptree_size) {

    // we now go through a(n over) complicated piece of (work) code to 
    // determine the amount of valid egress ifaces for these inputs. 
    int valid_ifaces = this->iface_num, i = 0;

    // if the local iface is not in the 'no forwarding' list, don't consider 
    // in the # of valid egress ifaces. why? the local iface is not an egress 
    // iface.
    std::set<int>::iterator it = this->no_forwarding.find((int) IFACE_LOCAL);
    if (it == this->no_forwarding.end()) {
        valid_ifaces -= 1;
        i = IFACE_LOCAL + 1;

        // std::cout << "RID_Router::get_valid_ifaces() [INFO] local link not included. valid_ifaces = " << valid_ifaces << std::endl;
    }

    // subtract the ifaces in the forwarding list, as well as those which 
    // don't have entries of size ptree_size
    for (; i < this->iface_num; i++) {

        it = this->no_forwarding.find(i);

        // remove ifaces in the 'no forwarding' list (i.e. not valid egress 
        // ifaces)
        if (it != this->no_forwarding.end()) {

            valid_ifaces--;

            // std::cout << "RID_Router::get_valid_ifaces() [INFO] link " << i << " included. valid_ifaces = " << valid_ifaces << std::endl;

        } else {

            // if ptree_iface.f_distribution[ptree_size] = 0.0, then it is not 
            // possible for a ptree to continue on ptree_iface
            __float080 f_dist = 0.0;
            if (ptree_size > 0)
                f_dist = this->fwd_table[i].f_distribution[ptree_size - 1];
            else
                f_dist = 1.0;

            // remove ifaces for which it's impossible to be in ptree of size 
            // ptree_size
            if (!(f_dist > 0.0)) {
                valid_ifaces--;

                // std::cout << "RID_Router::get_valid_ifaces() [INFO] LINK[" 
                //     << i << "].f_dist[" << (int) ptree_size << "] = 0.0 (" << f_dist 
                //     << "). valid_ifaces = " << valid_ifaces << std::endl;

                // for (int j = 0; j < this->f_max; j++)
                //     std::cout << "RID_Router::get_valid_ifaces() [INFO] F_DIST[" 
                //         << i << "].[" << j << "] = " << this->fwd_table[i].f_distribution[j] << std::endl;
            }
        }
    }

    // as a safety check, be sure valid_ifaces is at least 1.0
    if (valid_ifaces < 1) {
        valid_ifaces = 1;
    }

    return valid_ifaces;
}

__float080 RID_Router::calc_joint_log_prob(
    uint8_t ptree_size,
    uint8_t ptree_iface, 
    int * iface_pivots) {

    int valid_ifaces = get_valid_ifaces(ptree_size);

    // scale the calculated probabilities by the probability of having a 
    // request - of any size - enter this router
    __float080 log_prob = log(this->ingress_prob);
    log_prob += log(1.0 / ((__float080) valid_ifaces));

    for (int iface = IFACE_LOCAL; iface < this->iface_num; iface++) {

        // if iface is in the ptree of |F| = ptree_size, fetch the lpm_prob of 
        // ptree_size
        if (iface == ptree_iface) {

            log_prob += log(this->get_lpm_prob(iface, ptree_size, iface_pivots[iface]));

        } else {

            // otherwise, _iface is not in a tree (a FP may be generated 
            // at this router for the first time), so we should use 
            // ptree_size = 0
            log_prob += log(this->get_lpm_prob(iface, 0, iface_pivots[iface]));
        }
    }

    return log_prob;
}

int pivot_to_int(int * pivots, int pivot_size, long int pivot_init) {

    long int pivot = pivot_init;
    int p = 1, pp = p;

    // std::cout << "pivot_to_int() : [INFO] initial pivot : " << pivot << std::endl;

    for ( ; p < (pivot_size + 1); p++) {

        if (pivots[p - 1] > 9)
            pp++;

        pivot += ((long int) pivots[p - 1]) * ((long int) pow(10, pp++));

        // std::cout << "pivot_to_int() : [INFO] added pivot " << pivot << " from array ";
        // for (int i = 0; i < pivot_size; i++)
        //     std::cout << "[" << pivots[i] << "]";
        // std::cout << std::endl;
    }

    return pivot;
}

__float080 RID_Router::calc_cumulative_prob(
    __float080 * joint_prob_matrix,
    uint8_t iface, 
    uint8_t f,
    uint8_t mode,
    long int & pivot,
    __float080 & prob_new_events) {

    // printf("RID_Router::calc_cumulative_prob() : [IFACE : %d][|F| : %d][MODE : %s]\n", 
    //     iface, f, CUMULATIVE_PROB_MODE_STR[mode]);

    __float080 _prob = 0.0;
    __float080 _cumulative_prob = 0.0;

    // hold the pivots of the different iface arrays 
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int)); 

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;

    // FIXME: this is sort of hacky, but, well... it works ?
    uint8_t going_up = 0x00;

    // we fix the column (or row?) f for iface, and don't change it
    iface_pivots[iface] = f;

    // 
    int ptree_size = pivot;
    prob_new_events = 0.0;

    // the mode argument specifies the upper limit for |F| during the cumulative 
    // probability: 
    //  * if set to MODE_LI, we calculate P(L_{j,p} <= f_max), for 
    //    all j != iface. 
    //    why? no matter how large the matches for other ifaces are, the 
    //    local iface will always be chosen.
    //
    //  * if set to MODE_EI_EXCLUSIVE, we calculate P(L_{j,p} < f). 
    //
    //  * if set to MODE_EI_INCLUSIVE, we calculate P(L_{j,p} <= f). 
    //    why? because (...)
    //
    int f_limit = f;

    if (mode == MODE_EI_INCLUSIVE)
        f_limit = f + 1;

    if (mode == MODE_LI) 
        f_limit = this->f_max + 1;

    // if the mode is MODE_MIS, we create the iface_pivots in a special 
    // way : we keep iface_pivots[iface] = f, but only gather the joint prob 
    // values for iface_pivots with all indexes set to the same value
    if (mode == MODE_MIS) {

        if (this->iface_num < 3)
            return 0.0;

        for (int _f = 1; _f < (this->f_max + 1); _f++) {

            for (int i = 0; i < this->iface_num; i++) {

                iface_pivots[i] = _f;
            }

            // the pivots of the ifaces over which you CAN'T forward should be 
            // set to 0 (prob of that is 1.0)
            for (std::set<int>::iterator it = this->no_forwarding.begin();
                it != this->no_forwarding.end();
                ++it) {

                iface_pivots[(*it)] = 0;
            }

            iface_pivots[iface] = f;

            // printf("RID_Router::calc_cumulative_prob() : [MODE_MIS] joint_prob_matrix");        

            // for (int i = 0; i < this->iface_num; i++) 
            //     printf("[%d]", iface_pivots[i]);

            _prob = this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);
            _cumulative_prob += this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);

            // printf(" = %-.5LE (%-.5LE)\n", _prob, _cumulative_prob);
        }

        return _cumulative_prob;
    }

    while (curr_i > -1) {

        if (curr_i == (this->iface_num - 1)) {

            if (curr_i != iface) {

                // here we cycle through the the iface_pivots, only changing 
                // the pivot of the last iface. e.g. say we have 3 ifaces, and 
                // we're evaluating P(L_{0,p} < 5). one possible case for this 
                // cycle could be: 
                //  * -> [5][0][0] -> [5][0][1] -> ... -> [5][0][f - 1]
                // FIXME: this is where the MODE_LI mode can take effect
                for (int _f = 0; _f < f_limit; _f++) {

                    // update the |F| size pivot for iface i 
                    iface_pivots[this->iface_num - 1] = _f;
                    _prob = this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);

                    // convert the pivot array to an int and add it to a set. 
                    // this will keep track of events which are covered by 
                    // cumulative probabilities, and avoid the duplicate sum 
                    // of probability events
                    // FIXME: this could be a source of trouble
                    if (mode == MODE_EI_INCLUSIVE && _prob > 0.0) {

                        pivot = pivot_to_int(iface_pivots, this->iface_num, ptree_size);

                        std::set<long int>::iterator it = this->added_pivots.find(pivot);
                        if (it != this->added_pivots.end()) {

                            std::cout << "RID_Router::calc_cumulative_prob() : [INFO] event " 
                                << pivot << " has been considered before." << std::endl;

                        } else {

                            this->added_pivots.insert(pivot);

                            std::cout << "RID_Router::calc_cumulative_prob() : [INFO] added event " 
                                << pivot << "." << std::endl;

                            prob_new_events += _prob;
                        }
                    }

                    _cumulative_prob += _prob;

                    if (_prob > 0.0) {

                        printf("RID_Router::calc_cumulative_prob() : [%s] (1) joint_prob_matrix", CUMULATIVE_PROB_MODE_STR[mode]);        

                        for (int i = 0; i < this->iface_num; i++) 
                            printf("[%d]", iface_pivots[i]);

                        printf(" = %-.5LE (%-.5LE)\n", _prob, _cumulative_prob);
                    }
                }

            } else {

                _prob = this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);

                if (mode == MODE_EI_INCLUSIVE && _prob > 0.0) {

                    pivot = pivot_to_int(iface_pivots, this->iface_num, ptree_size);

                    std::set<long int>::iterator it = this->added_pivots.find(pivot);
                    if (it != this->added_pivots.end()) {

                        std::cout << "RID_Router::calc_cumulative_prob() : [INFO] event " 
                            << pivot << " has been considered before." << std::endl;

                    } else {

                        this->added_pivots.insert(pivot);

                        std::cout << "RID_Router::calc_cumulative_prob() : [INFO] added event " 
                            << pivot << "." << std::endl;

                        prob_new_events += _prob;
                    }
                }

                _cumulative_prob += _prob;
                going_up = 0x00;

                if (_prob > 0.0) {

                    printf("RID_Router::calc_cumulative_prob() : [%s] (2) joint_prob_matrix", CUMULATIVE_PROB_MODE_STR[mode]);        

                    for (int i = 0; i < this->iface_num; i++) 
                        printf("[%d]", iface_pivots[i]);

                    printf(" = %-.5LE (%-.5LE)\n", _prob, _cumulative_prob);
                }
            }
        }

        // extract the current pivot value for the curr_i level
        if (curr_i == iface && going_up) {
            
            // printf("RID_Router::calc_cumulative_prob() : [%s] checkpoint 1\n", CUMULATIVE_PROB_MODE_STR[mode]);
            curr_f = 0;

        } else {

            // printf("RID_Router::calc_cumulative_prob() : [%s] checkpoint 2\n", CUMULATIVE_PROB_MODE_STR[mode]);
            curr_f = iface_pivots[curr_i];
        }

        if (++curr_f < f_limit) {

            // printf("RID_Router::calc_cumulative_prob() : [%s] checkpoint 3 (curr_i = %d, curr_f = %d)\n", 
            //     CUMULATIVE_PROB_MODE_STR[mode], curr_i, curr_f);

            // FIXME: if we're only looking at EI events, we only need to consider 
            // those for which iface_pivots[IFACE_LOCAL] = 0. as soon as we 
            // increment iface_pivots[IFACE_LOCAL], get off the cycle.
            if ((mode == MODE_EI_EXCLUSIVE || mode == MODE_EI_INCLUSIVE) && curr_i == IFACE_LOCAL)
                break;

            // we're at some level curr_i, and we need to increment 
            // it's pivot counter (to curr_f)
            // e.g. we go from [5][0][0] -> [5][1][0]
            if (curr_i != iface) {
                iface_pivots[curr_i] = curr_f;
            }

            // after updating iface_pivots[curr_i], go back down 
            // to the level below:
            // old curr_i      new curr_i
            //     v               V
            // [5][0][0] -> [5][1][0]
            // we can then cycle through [5][1][0] -> [5][1][1] -> ... -> [5][1][f - 1]
            curr_i++;
            going_up = 0xFF;

        } else {

            // printf("RID_Router::calc_cumulative_prob() : [%s] checkpoint 4 (curr_i = %d, curr_f = %d)\n", 
            //     CUMULATIVE_PROB_MODE_STR[mode], curr_i, curr_f);

            // we've completed the iterations for the curr_i level. we 
            // now: 
            //  * reset iface_pivots[curr_i] to 0
            if (curr_i != iface) {

                // hack: the next time we do ++curr_f, curr_f will be 0
                curr_f = -1;
                iface_pivots[curr_i] = curr_f;

            }

            // order the function go up one level (i.e. to curr_i - 1) to increment 
            // iface_pivots[curr_i - 1]. 
            // e.g. we would go from [5][0][-1] -> [5][1][-1]
            // this only happens in the next cycle though
            if (--curr_i == IFACE_LOCAL && mode == MODE_LI)
                break;

            going_up = 0x00;
        }
    }

    free(iface_pivots);

    return _cumulative_prob;
}

int RID_Router::calc_iface_events_pmf(__float080 * joint_prob_matrix, int ptree_size) {

    long int pivot = 0;
    __float080 prob_new_events = 0.0;

    // LOCAL MATCH : call calc_cumulative_prob in MODE_LI mode
    for (uint8_t _f = this->f_max; _f > 0; _f--) {

        // printf("RID_Router::calc_iface_events_pmf() : calc_cumulative_prob(..., %d, %d, %s)\n", 
        //     IFACE_LOCAL, _f, CUMULATIVE_PROB_MODE_STR[MODE_LI]);
        this->iface_events_pmf[EVENT_LLM] += this->calc_cumulative_prob(joint_prob_matrix, IFACE_LOCAL, _f, MODE_LI, pivot, prob_new_events);
    }

    std::cout << "RID_Router::calc_iface_events_pmf() : [INFO]"
        << "\n\t P('EVENT_LLM' [ptree_size = " << (int) ptree_size << "]') = " << this->iface_events_pmf[EVENT_LLM]
        << std::endl;

    // NO MATCHES : just call get_joint_lpm_prob() with iface_pivots = {0, 0, ..., 0}
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    this->iface_events_pmf[EVENT_NLM] += this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);

    // MULTIPLE & SINGLE IFACE MATCHES 
    __float080 prob_single = 0.0;
    __float080 prob_multiple = 0.0;

    for (uint8_t iface = 1; iface < this->iface_num; iface++) {
        for (uint8_t f = this->f_max; f > 0; f--) {

            // reset the pivot 'value-result' arg for calc_cumulative_prob()
            pivot = ptree_size;

            // probability of choosing iface being *over others*, with an entry 
            // of size f.
            prob_single = this->calc_cumulative_prob(joint_prob_matrix, iface, f, MODE_EI_EXCLUSIVE, pivot, prob_new_events);
            // in some cases, we may want the probability of iface being chosen, 
            // even if together with any other iface (e.g. if a flooding strategy 
            // is in place). 
            prob_multiple = this->calc_cumulative_prob(joint_prob_matrix, iface, f, MODE_EI_INCLUSIVE, pivot, prob_new_events);

            // keep track of the total EI probability as well (sanity check)
            if (prob_single > 0.0) {

                std::cout << "RID_Router::calc_iface_events_pmf() : [INFO]"
                    << "\n\t P('EGRESS_IFACE_EXCLUSIVE' [" << (int) iface << "][" << (int) f << "]') = " << prob_single
                    << std::endl;
            }

            if (prob_multiple > 0.0) {

                std::cout << "RID_Router::calc_iface_events_pmf() : [INFO]"
                    << "\n\t P('EGRESS_IFACE_INCLUSIVE' [" << (int) iface << "][" << (int) f << "]') = " << prob_multiple
                    << std::endl;
            }

            // depending on the way RID routers handle multiple matches, we 
            // set the probability of choosing the iface differently.
            if (this->mm_mode == MMH_FLOOD) {

                this->egress_iface_prob[iface][f] += prob_multiple;

            } else if (this->mm_mode == MMH_RANDOM) {

                // FIXME: i've identified a problem: if the outdegree 
                // is > 2, then events as [2][2][1] should be divided 
                // by 1/2, while events such as [2][2][2] should be divided 
                // by 1/3. this only fails for MMH_RANDOM.
                int valid_ifaces = get_valid_ifaces(ptree_size);
                this->egress_iface_prob[iface][f] += (prob_single + ((1.0 / (__float080) valid_ifaces) * (prob_multiple - prob_single)));


                if ((prob_multiple - prob_single) > 0.0) {

                    std::cout << "RID_Router::calc_iface_events_pmf() : [INFO] "
                        << "\n\t [MMH_RANDOM] P('EGRESS_IFACE[" << (int) iface << "][" << (int) f << "][ptree_size = " << ptree_size << "]') += " 
                        << (prob_single + ((1.0 / (__float080) valid_ifaces) * (prob_multiple - prob_single)))
                        << std::endl;

                    std::cout << "RID_Router::calc_iface_events_pmf() : [INFO] "
                        << "\n\t [MMH_RANDOM] P('EGRESS_IFACE[" << (int) iface << "][" << (int) f << "][ptree_size = " << ptree_size << "]') = " << this->egress_iface_prob[iface][f]
                        << std::endl;

                    std::cout << "RID_Router::calc_iface_events_pmf() : [INFO] "
                        << "\n\t valid_ifaces = " << valid_ifaces << std::endl;
                }

            } else {

                this->egress_iface_prob[iface][f] += prob_single;
            }

            // regardless of the way multiple matches are handled, we track 
            // the probability of an exclusive iface choice event correctly
            this->iface_events_pmf[EVENT_SLM] += prob_single;

            if (prob_new_events > 0.0) {

                this->iface_events_pmf[EVENT_MLM] += (prob_new_events - prob_single);
                std::cout << "RID_Router::calc_iface_events_pmf() : [INFO] P('EVENT_MLM' [ptree_size = " << ptree_size << "]) = " 
                    << this->iface_events_pmf[EVENT_MLM] << std::endl;
            }

            // keep track of the total EI probability as well (sanity check)
            if (prob_single > 0.0 || prob_multiple > 0.0) {

                std::cout << "RID_Router::calc_iface_events_pmf() : [INFO] "
                    << "\n\t P('EGRESS_IFACE[" << (int) iface << "][" << (int) f << "]') = " << this->egress_iface_prob[iface][f]
                    << std::endl;
            }
        }
    }

    // clean up allocated memory
    free(iface_pivots);

    return 0;
}

__float080 RID_Router::get_iface_events_prob(uint8_t event) {

    if (event < 0 || event > EVENT_SLM) {

        fprintf(stderr, "RID_Router::get_iface_events_prob() : unknown event nr. (%d)\n", event);

        return 0.0;
    }

    return this->iface_events_pmf[event];
}

__float080 * RID_Router::get_egress_ptree_prob(uint8_t iface) {

    return this->egress_ptree_prob[iface];
}

__float080 RID_Router::get_egress_iface_prob(uint8_t iface) {

    __float080 iface_prob = 0.0;

    for (int f = 0; f < (this->f_max + 1); f++)
        iface_prob += this->egress_iface_prob[iface][f];

    return iface_prob;
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
