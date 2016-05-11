#include "rid-router.h"

const char * EVENT_STR[] = { "NIS", "MIS", "LI", "EI"};
const char * CUMULATIVE_PROB_MODE_STR[] = { "MODE_EI", "MODE_MIS", "MODE_LI"};

RID_Router::RID_Router(
    uint8_t access_tree_index, 
    uint8_t height, 
    uint8_t width,
    uint32_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max,
    uint16_t bf_size) {

    this->access_tree_index = access_tree_index;
    this->height = height;
    this->width = width;
    this->starting_router = false;

    this->total_joint_prob = 0.0;

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

    // initialize the ingress_size_pmf[] array
    this->ingress_size_pmf = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));

    // initialize the egress_size_pmf[iface][] array
    this->egress_size_pmf = (__float080 **) calloc(this->iface_num, sizeof(__float080));
    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++)
        this->egress_size_pmf[_iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
}

RID_Router::~RID_Router() {

    // free forwarding table
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
        free(this->egress_size_pmf[_iface]);
}

int RID_Router::forward(
    uint8_t request_size,     
    uint8_t ingress_iface,              
    int * tp_sizes,
    __float080 * ingress_probs,                 
    __float080 * f_r_distribution) {

    this->iface_ingress = ingress_iface;

    for (uint8_t _ptree_size = 0; _ptree_size < this->f_max + 1; _ptree_size++) {
        
        this->ingress_size_pmf[_ptree_size] = ingress_probs[_ptree_size];

        // printf("RID_Router::forward() : ingress_probs[%d] = %-.5LE\n", 
        //     _ptree_size, (__float080) this->ingress_size_pmf[_ptree_size]);
    }

    // printf("RID_Router::forward() : 1 / (iface_num = %d) = %-.5LE\n", 
    //     this->iface_num - 1, 
    //     (__float080) 1.0 / (__float080) (this->iface_num - 1));
    
    // 1) calculate the distribution for the L_{i,p} random variables : for 
    // each iface i and 'prefix tree' size p
    if (calc_lpm_pmf(
        request_size, 
        ingress_iface, 
        tp_sizes,              // P('prefix tree size = p')
        f_r_distribution) < 0) {

        return -1;
    }

    for (uint8_t _iface = IFACE_LOCAL; _iface < this->iface_num; _iface++)
        this->print_lpm_pmf(_iface);

    // 2) calculate the joint probability distribution of the 
    // random variables L_{0,p} x L_{1,p} x ... x L_{|R|_{max},p}.

    // 2.1) initialize a single array to compute the joint prob distribution for 
    // all 'prefix tree' sizes ptree_size
    __float080 * _joint_lpm_matrix = NULL;
    init_joint_lpm_pmf(&(_joint_lpm_matrix));

    // 2.2) compute the joint distr. for p, considering each diff. iface 
    // as the iface associated with the prefix tree 
    for (uint8_t _ptree_size = 0; _ptree_size < (this->f_max + 1); _ptree_size++) {

        // clear the array for each diff. value of p
        clear_joint_lpm_pmf(&(_joint_lpm_matrix));

        for (uint8_t _ptree_iface = 0; _ptree_iface < this->iface_num; _ptree_iface++) {

            if (_ptree_iface == this->iface_ingress)
                continue;

            // printf("RID_Router::forward() : [ptree_size %d][iface %d]\n", 
            //     _ptree_size, _ptree_iface);

            calc_joint_lpm_pmf(_joint_lpm_matrix, _ptree_size, _ptree_iface);
        }

        // printf("RID_Router::forward() : total_joint_prob = %-.5LE\n", 
        //     (__float080) this->total_joint_prob);

        // 3) calc the prob distribution for ifaces, for each ptree_size
        if (calc_iface_events_pmf(_joint_lpm_matrix) < 0) {

            return -1;
        }
    }

    // FIXME : i don't know if this is the correct way of doing this...
    // if the request is going up, we relay NIS cases to the router upstream. 
    // once we reach the top of the topology, we give up relaying and use 
    // a 'fallback' to origin
    if (ingress_iface != IFACE_UPSTREAM 
        && (this->height > 0 && this->width > 0)) {
        this->egress_size_pmf[IFACE_UPSTREAM][0] = this->iface_events_pmf[EVENT_NIS];
    }

    this->print_iface_events_pmf();
    this->print_egress_size_pmf();

    // think it's safe to do this? think again...
    free(_joint_lpm_matrix);

    return 0;
}

int RID_Router::add_fwd_table_entry(
    int iface, 
    __float080 iface_proportion, 
    __float080 * f_distribution) {

    if (iface >= IFACE_LOCAL && iface < iface_num) {

        // this marks the entry as valid, as it becomes != -1
        this->fwd_table[iface].iface = iface;
        this->fwd_table[iface].iface_proportion = iface_proportion;
        this->fwd_table[iface].f_distribution = f_distribution;
        this->fwd_table[iface].next_hop = NULL;

        // printf("RID_Router::add_fwd_table_entry() : added fwd entry to r[%d][%d]:"\
        //     "\n\t[iface] = %d"\
        //     "\n\t[iface_proportion] = %-.5LE\n", 
        //     this->height, this->width, this->fwd_table[iface].iface, this->fwd_table[iface].iface_proportion = iface_proportion);

    } else {

        fprintf(
            stderr, 
            "RID_Router::add_fwd_table_entry() : invalid iface index (%d, "\
                "should be in [0, 1, ..., %d])\n", 
            iface, (this->iface_num - 1));
    }

    return 0;
}

void RID_Router::set_fwd_table_next_hop(uint8_t iface, RID_Router * next_hop_router) {
    this->fwd_table[iface].next_hop = next_hop_router;
}

RID_Router * RID_Router::get_fwd_table_next_hop(uint8_t iface) {
    return this->fwd_table[iface].next_hop;
}

void RID_Router::set_fwd_table(RID_Router::fwd_table_row * fwd_table) { 
    this->fwd_table = fwd_table; 
}

RID_Router::fwd_table_row * RID_Router::get_fwd_table() {

    return this->fwd_table;
}

uint32_t RID_Router::get_fwd_table_size() { 
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

void RID_Router::set_as_starting_router() {

    this->starting_router = true;
}

bool RID_Router::is_starting_router() {

    return this->starting_router;
}

__float080 fp_rate(__float080 m, __float080 n, __float080 k, __float080 c) {

    printf("fp_rate() :"\
        "\n\tm = %-.5LE"\
        "\n\tn = %-.5LE"\
        "\n\tk = %-.5LE"\
        "\n\t|F\\R| = %-.5LE"\
        "\n\tc = %-.5LE\n",
        m, n, k, c, (__float080) pow((1.0 - exp(-((n / m) * k))), k * c));


    return pow((1.0 - exp(-((n / m) * k))), k * c);
}

int RID_Router::get_log_fp_rates(
    __float080 m, 
    __float080 n, 
    __float080 k,
    __float080 iface_proportion,
    __float080 * f_distribution, 
    __float080 * f_r_distribution,
    __float080 * log_fp_rates) {        // function fills this array

    __float080 n_entries = 0.0;
    __float080 f_entries = 0.0; 

    for (int f = 0; f < this->f_max; f++) {

        f_entries = (__float080) this->fwd_table_size * iface_proportion * f_distribution[f];

        for (int i = 0; i <= f; i++) {

            n_entries = (__float080) f_entries * f_r_distribution[i];
            log_fp_rates[f] += n_entries * log(1.0 - fp_rate(m, n, k, (__float080) (i + 1)));
        }
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
 *          <iface, prefix tree size> pairs in the router
 *
 * the random variable L_{i,p} represents the longest match reported by an 
 * iface i, for a prefix tree of size p. P(L_{i,p} = f) is the probability 
 * of getting a longest match equal to f.
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
 * positive matches can happen due to (1) false positive matches and (2) true 
 * positive matches. in their turn, FP matches can happen in 2 ways:
 *  -# 'local' FP matches, triggered for the first time at this router 
 *      (in other words, a request found a new tree)
 *  -# FP matches coming from a previous router, due to a 'prefix tree binding'
 * 
 * in addition to FP matches, P(L_{i,p} = f) is influenced by TP matches, provided 
 * as input to the function. unlike FPs, true positive info is 
 * deterministic (as opposed to "probabilistic").
 * 
 * \param   fp_size_prob    probabilities of the FP match sizes that led the 
 *                          request to this router.
 *
 * \return  1 if bit i is set, 0 if not set.
 */
int RID_Router::calc_lpm_pmf(
    uint8_t request_size, 
    uint8_t ingress_iface,
    int * tp_sizes, 
    __float080 * f_r_distribution) {

    // probability of having an iface as part of a 'prefix tree' of size |F| = f

    // speed things up by keeping the (~FP5 * ~FP4 * ~FP4 * ... * ~FPi)
    // probabilities in an array
    __float080 * _log_prob_not_fp = (__float080 *) calloc(request_size + 1, sizeof(__float080));
    __float080 * _log_fp_rates = (__float080 *) calloc(this->f_max, sizeof(__float080));

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
            
            // let's assume that, according to LPM rules, we don't forward 
            // packets over the ingress iface
            if (_iface == ingress_iface) {

                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[0] = 1.0;
            }

            // FIXME: as a special fixme, if the next hop is 'NULL', treat it 
            // in the same way as the ingress iface
            if (_iface > IFACE_LOCAL && this->get_fwd_table_next_hop(_iface) == NULL) {

                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[0] = 1.0;
            }
        }

        // thus, skip everything else if _iface is the ingress iface
        if (_iface == ingress_iface) 
            continue;

        if (_iface > IFACE_LOCAL && this->get_fwd_table_next_hop(_iface) == NULL) 
            continue;

        // get the 1 x |R|_max matrix with the probs of *NOT* having FPs among  
        // entries of size f in _iface. these probs are given in log() 
        // form, since we'll be dealing with small numbers.
        __float080 _k = (log(2) * ((__float080) this->bf_size)) / ((__float080) request_size);

        if (get_log_fp_rates(
            (__float080) (this->bf_size), (__float080) request_size, _k,   // BF parameters
            this->fwd_table[_iface].iface_proportion,   // % of entries for _iface
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

        // if (_iface == 0)
        //     printf("RID_Router::calc_lpm_pmf() : P( |FP|(_iface = %d) <= %d) = %-.5LE\n", _iface, this->f_max, (__float080) exp(_log_prob_not_fp[this->f_max]));

        for (int _f = (this->f_max - 1); _f >= 0; _f--) {
            _log_prob_not_fp[_f] = _log_prob_not_fp[_f + 1] + _log_fp_rates[_f];

            // if (_iface == 0) {

            //     printf("RID_Router::calc_lpm_pmf() : P(~|FP|(_iface = %d) == %d) = %-.5LE\n", _iface, _f + 1, (__float080) exp(_log_fp_rates[_f]));
            //     printf("RID_Router::calc_lpm_pmf() : P( |FP|(_iface = %d) == %d) = %-.5LE\n", _iface, _f + 1, (__float080) 1.0 - (__float080) exp(_log_prob_not_fp[_f]));
            //     printf("RID_Router::calc_lpm_pmf() : P( |FP|(_iface = %d) <= %d) = %-.5LE\n", _iface, _f, (__float080) exp(_log_prob_not_fp[_f]));
            // }
        }

        // start with P(L_{_iface,_ptree_size} = |R|_max AND 
        // _ptree_size = |R|_max), go down from there...
        for (int _ptree_size = this->f_max; _ptree_size >= 0; _ptree_size--) {

            for (int _f = this->f_max; _f >= 0; _f--) {

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
                                (__float080) ((this->fwd_table[_iface].iface_proportion > 0.0) ? 0.0 : 1.0);

                        } else {

                            // if _f >= _ptree_size, the longest match will *AT 
                            // LEAST* be _f. the expression becomes: 
                            // (~FP5 * ... * ~FP[_f + 1])   : no FPs with |FP| > _f
                            // * ((1 - ~FP[_f]) + ~FP[_f])  : it doesn't matter if 
                            //                      a 'local' FP match with 
                            //                      |FP| > _f happens: we 
                            //                      will have a match for size _f 
                            //                      anyway
                            this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                exp(_log_prob_not_fp[_f]);
                        }

                    } else {

                        // 2) no TP exists for |F| = _f

                        // so, the expression becomes: 
                        // (~FP5 * ... * ~FP[_f + 1])   : no FPs with |FP| > _f
                        // * (1 - ~FP[_f])              : a 'local' FP match with 
                        //                          |FP| > _f *MUST* happen if the 
                        //                          iface is *NOT* on a ptree
                        if (_f > 0) {

                            if (_f < _ptree_size) {

                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 0.0;

                                // if (_iface == 0) {

                                //     printf("RID_Router::calc_lpm_pmf() :|F| < p, P(%d == %d) = %-.5LE\n", 
                                //         _f, _ptree_size, (__float080) 0.0);
                                // }

                            } else if (_f == _ptree_size) {


                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    (__float080) exp(_log_prob_not_fp[_f]) * ((this->fwd_table[_iface].iface_proportion > 0.0) ? 1.0 : 0.0);

                                // if (_iface == 0) {

                                //     printf("RID_Router::calc_lpm_pmf() :|F| == p, P(%d == %d) = %-.5LE\n", 
                                //         _f, 
                                //         _ptree_size, 
                                //         (__float080) exp(_log_prob_not_fp[_f]) * (__float080) ((this->fwd_table[_iface].iface_proportion > 0.0) ? (__float080) 1.0 : (__float080) 0.0));

                                //     printf("RID_Router::calc_lpm_pmf() : P(|F|_{%d} = %d) = %-.5LE, (ternary op.) = %-.5LE\n", 
                                //         _iface,
                                //         _f,
                                //         (__float080) this->fwd_table[_iface].iface_proportion,
                                //         (__float080) ((this->fwd_table[_iface].iface_proportion > 0.0) ? 1.0 : 0.0));
                                // }

                            } else {

                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    exp(_log_prob_not_fp[_f]) * (1.0 - exp(_log_fp_rates[_f - 1]));

                                // if (_iface == 0) {

                                //     printf("RID_Router::calc_lpm_pmf() :|F| > p, P(%d == %d) = %-.5LE\n", 
                                //         _f, _ptree_size, (__float080) exp(_log_prob_not_fp[_f]) * (1.0 - exp(_log_fp_rates[_f - 1])));
                                // }
                            }

                        } else {

                            if (_f < _ptree_size) {

                                // FIXME: this seems hack-ish : here i check if 
                                // _iface is associated with any entries. if it 
                                // is, then
                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    (__float080) ((this->fwd_table[_iface].iface_proportion > 0.0) ? 0.0 : 1.0);

                            } else {

                                this->lpm_pmf[_iface][_ptree_size].lpm_pmf_prob[_f] = 
                                    exp(_log_prob_not_fp[_f]);
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
 * \return  0 if everything goes fine, an error code < 0 if errors occur.
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
                _log_prob = exp(_log_prob) * this->ingress_size_pmf[ptree_size];

                // set this value in the joint_lpm_size_matrix 
                // note that this function adds exp(log_prob) to the value 
                // already held by joint_lpm_size_matrix in the position indexed 
                // by iface_pivots
                this->add_joint_lpm_prob(joint_prob_matrix, iface_pivots, _log_prob);
                this->total_joint_prob += _log_prob;

                // int _sum = 0;
                // for (int _iface = 0; _iface < this->iface_num; _iface++) 
                //     _sum += iface_pivots[_iface];

                // if (_sum == 0) {

                //     printf("RID_Router::calc_joint_lpm_pmf() : [ptree_size = %d][ptree_iface = %d] joint_prob_matrix",
                //         ptree_size, ptree_iface);        

                //     for (int _iface = 0; _iface < this->iface_num; _iface++) 
                //         printf("[%d]", iface_pivots[_iface]);

                //     printf(" = %-.5LE (%-.5LE)\n", _log_prob, this->total_joint_prob);
                // }
            }
        }

        // extract the current pivot value for the curr_i level
        curr_f = iface_pivots[curr_i];

        if (++curr_f < (this->f_max + 1)) {

            // we're at some level curr_i, and we need to increment 
            // it's pivot counter (to curr_f)
            iface_pivots[curr_i] = curr_f;

            // after updating iface_pivots[curr_i], go back down 
            // to the level below.
            curr_i++;

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

__float080 RID_Router::calc_joint_log_prob(
    uint8_t ptree_size,
    uint8_t ptree_iface, 
    int * iface_pivots) {

    __float080 _log_prob = 0.0;

    // // all probabilities are scaled by (1 / iface_num). why '-1'? we don't 
    // // count with the ingress iface
    if (this->starting_router == false)
        _log_prob = log((1.0) / (__float080) (this->iface_num - 1));
    else
        _log_prob = log((1.0) / (__float080) (this->iface_num));
    // _log_prob = log((1.0) / (__float080) (this->iface_num));

    for (uint8_t _iface = 0; _iface < this->iface_num; _iface++) {

        // if _iface is in the ptree of ptree_size, fetch the lpm_prob of 
        // ptree_size
        if (_iface == ptree_iface)
            _log_prob += log(this->get_lpm_prob(_iface, ptree_size, iface_pivots[_iface]));
        else
            // otherwise, _iface is not in a tree, so we should use 
            // ptree_size = 0
            _log_prob += log(this->get_lpm_prob(_iface, 0, iface_pivots[_iface]));
    }

    return _log_prob;
}

__float080 RID_Router::calc_cumulative_prob(
    __float080 * joint_prob_matrix,
    uint8_t iface, 
    uint8_t f,
    uint8_t mode) {

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

    // the mode argument specifies the upper limit for the cumulative 
    // probability: 
    //  * if set to MODE_LI, we calculate P(L_{j,p} <= f_max), for 
    //      all j =/= iface. 
    //  * if set to MODE_EI, we calculate P(L_{j,p} < f). 
    //
    // the MODE_LI mode is usually set when calculating the probability 
    // of the LI (LOCAL) interface event
    int _upper_limit = f;

    if (mode == MODE_LI) 
        _upper_limit = this->f_max + 1;

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

            iface_pivots[this->iface_ingress] = 0;
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
                for (int _f = 0; _f < _upper_limit; _f++) {

                    // update the |F| size pivot for iface i 
                    iface_pivots[this->iface_num - 1] = _f;
                    _prob = this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);
                    _cumulative_prob += _prob;

                    // printf("RID_Router::calc_cumulative_prob() : [%s] (1) joint_prob_matrix", CUMULATIVE_PROB_MODE_STR[mode]);        

                    // for (int i = 0; i < this->iface_num; i++) 
                    //     printf("[%d]", iface_pivots[i]);

                    // printf(" = %-.5LE (%-.5LE)\n", _prob, _cumulative_prob);
                }

            } else {

                _prob = this->get_joint_lpm_prob(joint_prob_matrix, iface_pivots);
                _cumulative_prob += _prob;
                going_up = 0x00;

                // printf("RID_Router::calc_cumulative_prob() : [%s] (2) joint_prob_matrix", CUMULATIVE_PROB_MODE_STR[mode]);        

                // for (int i = 0; i < this->iface_num; i++) 
                //     printf("[%d]", iface_pivots[i]);

                // printf(" = %-.5LE (%-.5LE)\n", _prob, _cumulative_prob);
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

        if (++curr_f < _upper_limit) {

            // printf("RID_Router::calc_cumulative_prob() : [%s] checkpoint 3 (curr_i = %d, curr_f = %d)\n", 
            //     CUMULATIVE_PROB_MODE_STR[mode], curr_i, curr_f);

            // FIXME: if we're only looking at EI events, we only need to consider 
            // those for which iface_pivots[IFACE_LOCAL] = 0. as soon as we 
            // increment iface_pivots[IFACE_LOCAL], get off the cycle.
            if (mode == MODE_EI && curr_i == IFACE_LOCAL)
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

int RID_Router::calc_iface_events_pmf(
    __float080 * joint_prob_matrix) {

    __float080 _prob = 0.0;

    // LI event : call calc_cumulative_prob in MODE_LI mode
    for (uint8_t _f = this->f_max; _f > 0; _f--) {

        // printf("RID_Router::calc_iface_events_pmf() : calc_cumulative_prob(..., %d, %d, %s)\n", 
        //     IFACE_LOCAL, _f, CUMULATIVE_PROB_MODE_STR[MODE_LI]);
        this->iface_events_pmf[EVENT_LI] += this->calc_cumulative_prob(joint_prob_matrix, IFACE_LOCAL, _f, MODE_LI);
    }

    // NIS event : just call get_joint_lpm_prob() with iface_pivots = {0, 0, ..., 0}
    int * _iface_pivots = (int *) calloc(this->iface_num, sizeof(int));
    this->iface_events_pmf[EVENT_NIS] += this->get_joint_lpm_prob(joint_prob_matrix, _iface_pivots);

    // MIS events : call calc_cumulative_prob in MODE_MIS mode, 
    // setting IFACE_LOCAL to f = 0
    if (this->iface_num < 3)
        this->iface_events_pmf[EVENT_MIS] = 0.0;
    else
        this->iface_events_pmf[EVENT_MIS] += this->calc_cumulative_prob(joint_prob_matrix, IFACE_LOCAL, 0, MODE_MIS);

    // EI events : this is a bit more complicated, as we need to fill 
    // the egress_size_pmf array, for each iface
    for (uint8_t _iface = 1; _iface < this->iface_num; _iface++) {

        for (uint8_t _f = this->f_max; _f > 0; _f--) {

            // P(I = i AND |F| = f) : this is what we want for the egress probs. 
            _prob = this->calc_cumulative_prob(joint_prob_matrix, _iface, _f, MODE_EI);
            this->egress_size_pmf[_iface][_f] += _prob;

            // keep track of the total EI probability as well (sanity check)
            this->iface_events_pmf[EVENT_EI] += _prob;
            // printf("RID_Router::calc_cumulative_prob() : egress_size_pmf[%d][%d] = %-.5LE (%-.5LE) (P(EVENT_EI) = %-.5LE)\n", 
            //     _iface, _f, _prob, this->egress_size_pmf[_iface][_f], this->iface_events_pmf[EVENT_EI]);
        }
    }

    // clean up allocated memory
    free(_iface_pivots);

    return 0;
}

__float080 RID_Router::get_iface_events_prob(uint8_t event) {

    if (event < 0 || event > EVENT_EI) {

        fprintf(stderr, "RID_Router::get_iface_events_prob() : unknown event nr. (%d)\n", event);

        return 0.0;
    }

    return this->iface_events_pmf[event];
}

__float080 RID_Router::get_egress_size_prob(uint8_t iface, uint8_t f) {

    return this->egress_size_pmf[iface][f];
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
    for (uint8_t _event = 0; _event < EVENT_NUM + 1; _event++)
        printf("------------");

    printf("\n[EVENT] :  |");

    for (int _event = 0; _event < EVENT_NUM; _event++)
        printf(" %-11s|", EVENT_STR[_event]);

    printf("  SUM P(e)  |");

    printf("\n-----------------");
    for (uint8_t _event = 0; _event < EVENT_NUM + 1; _event++)
        printf("------------");

    __float080 _prob = 0.0;
    __float080 _prob_total = 0.0;

    printf("\n[P(e)] :   |");

    for (uint8_t _event = 0; _event < EVENT_NUM; _event++) {

        _prob = this->iface_events_pmf[_event];
        _prob_total += _prob;

        printf(" %-.5LE|", _prob);
    }

    printf(" %-.5LE|", _prob_total);

    printf("\n-----------------");
    for (uint8_t _event = 0; _event < EVENT_NUM + 1; _event++)
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
void RID_Router::print_egress_size_pmf() {

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

            _prob = this->egress_size_pmf[_iface][_f];
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
