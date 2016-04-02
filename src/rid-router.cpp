#include "rid-router.h"

RID_Router::RID_Router(
    uint8_t access_tree_index, 
    uint8_t height, 
    uint8_t width,
    uint32_t fwd_table_size,
    uint8_t iface_num,
    uint8_t f_max) {

    this->access_tree_index = access_tree_index;
    this->height = height;
    this->width = width;

    this->fwd_table_size = fwd_table_size;
    this->iface_num = iface_num;
    this->f_max = f_max;

    // initialize forwarding table
    this->fwd_table = (RID_Router::fwd_table_row *) calloc(iface_num, sizeof(RID_Router::fwd_table_row)); 
    // to indicate an empty entry in the table, set iface to -1
    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++)
        this->fwd_table[i].iface = -1;
}

RID_Router::~RID_Router() {

    // forwarding table
    free(this->fwd_table);

    // pmf for F_i random variable
    for (int i = 0; i < this->iface_num; i++)
        free(this->f_pmf[i]);

    free(this->f_pmf);

    // joint probability distribution
    free(this->joint_f_pmf);

    // pmf for I_i random variable
    free(this->lpm_iface_pmf);
}

int RID_Router::forward(
    uint8_t request_size,     
    uint8_t ingress_iface,              
    int * tp_sizes,                 
    __float080 * f_r_distribution) {
    
    // 1) calculate the PMFs for the F_i random variables
    if (calc_f_pmf(request_size, ingress_iface, tp_sizes, f_r_distribution) < 0) {

        return -1;
    }

    this->print_f_pmf();

    // 2) calculate the joint probability distribution of the 
    // random variables {F_0, F_1, ..., F_n}.
    if (calc_joint_f_pmf() < 0) {

        return -1;
    }

    //this->get_path_state()->print_joint_f_pmf();

    // 3) finally, calculate the PMF for a random variable I, which 
    // represents the interface chosen by the LPM engine of the router, 
    // accounting for the presence of FPs. I can take the values {-1, 0, 
    // 1, ..., iface_num - 1}. '-1' is used for a 'no match' event.

    // the PMF for I is calculated from the the joint probability 
    // distribution of {F_0, F_1, ..., F_n}
    if (calc_lpm_iface_pmf() < 0) {

        return -1;
    }

    this->print_lpm_iface_pmf();

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

__float080 fp_rate(__float080 m, __float080 n, __float080 k, __float080 c) {

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

__float080 * RID_Router::init_f_pmf(uint8_t iface) {

    // f_pmf is a 2D array : it may not be initialized when init_f_pmf() is 
    // first called
    if (this->f_pmf == NULL)
        this->f_pmf = (__float080 **) calloc(this->iface_num, sizeof(__float080));

    this->f_pmf[iface] = (__float080 *) calloc(this->f_max + 1, sizeof(__float080));
    return this->f_pmf[iface];
}

void RID_Router::set_f_pmf(int iface, int f, __float080 prob) {
    this->f_pmf[iface][f] = prob;
}

__float080 * RID_Router::get_f_pmf(uint8_t iface) { 
    return this->f_pmf[iface]; 
}

__float080 RID_Router::get_f_pmf(uint8_t iface, uint8_t f) { 
    return this->f_pmf[iface][f]; 
}

void RID_Router::print_f_pmf() {

    printf("\n------------------");
    for (uint8_t i = 0; i < (this->f_max + 2); i++)
        printf("-------------");

    printf("\n[|F|]          : |");

    for (uint8_t i = 0; i < (this->f_max + 1); i++)
        printf(" %-11d|", i);

    printf(" SUM P(F_i) |");

    printf("\n------------------");
    for (uint8_t i = 0; i < (this->f_max + 2); i++)
        printf("-------------");

    __float080 subtotal = 0.0;
    __float080 f_pmf_prob = 0.0;

    for (uint8_t i = 0; i < (this->iface_num); i++) {

        printf("\n[P(F_%d = |F|)] : |", i);

        for (uint8_t f = 0; f < (this->f_max + 1); f++) {

            f_pmf_prob = this->get_f_pmf(i, f);
            printf(" %-.5LE|", f_pmf_prob);

            subtotal += this->get_f_pmf(i, f);
        }

        printf(" %-.5LE|", subtotal);

        // reset the subtotal accumulator
        subtotal = 0.0;
    }

    printf("\n------------------");
    for (uint8_t i = 0; i < (this->f_max + 2); i++)
        printf("-------------");

    printf("\n");
}

__float080 * RID_Router::init_joint_f_pmf() {

    // initialize joint_f_pmf 'matrix'. since n-dimensional arrays in 
    // C/C++ are initialized in row major order, we access the 
    // contents of joint_f_pmf by 'smartly' dereferencing the pointer (wooo, 
    // watch out, we've got a rocket scientist over here...),
    this->joint_f_pmf = (__float080 *) calloc(pow((this->f_max + 1), this->iface_num), sizeof(__float080));

    return this->joint_f_pmf;
}

void RID_Router::set_joint_f_pmf(
    int * iface_pivots,
    __float080 value) {

    uint32_t joint_f_pmf_add = 0;

    // arrays in C/C++ are stored in row major order. what does it 
    // mean? e.g. imagine that iface_num = 3, which yields a 3D matrix. 
    // to access element joint_f_pmf[1][3][5], we access the address 
    // this->joint_f_pmf + (1 * ((f_max + 1)^2)) + (3 * ((f_max + 1)^1)) + (5 * ((f_max + 1)^0))
    for (uint8_t i = 0; i < this->iface_num; i++) {

        joint_f_pmf_add += iface_pivots[i] * pow((this->f_max + 1), i);
    }

    *(this->joint_f_pmf + joint_f_pmf_add) = value;
}

__float080 RID_Router::get_joint_f_pmf(int * iface_pivots) {

    uint32_t joint_f_pmf_add = 0;

    // arrays in C/C++ are stored in row major order. what does it 
    // mean? e.g. imagine that iface_num = 3, which yields a 3D matrix. 
    // to access element joint_f_pmf[1][3][5], we access the address 
    // this->joint_f_pmf + (1 * ((f_max + 1)^2)) + (3 * ((f_max + 1)^1)) + (5 * ((f_max + 1)^0))
    for (int i = 0; i < this->iface_num; i++) {

        joint_f_pmf_add += iface_pivots[i] * pow((this->f_max + 1), i);
    }

    return *(this->joint_f_pmf + joint_f_pmf_add);
}

void RID_Router::print_joint_f_pmf() {

    if (this->iface_num > 3) {

        fprintf(stderr, "RID_Router::print_joint_f_pmf() : cannot print"\
            " 3D LPM PMF matrix, iface_num != 3 (%d)\n", 
            this->iface_num);

        return;
    }

    int iface_pivots[3] = {0, 0, 0};
    __float080 total = 0.0;

    for (int f_0 = 0; f_0 < (this->f_max + 1); f_0++) {

        iface_pivots[0] = f_0;

        // print f_max 2D matrices, one fore each possible value of F_0
        printf("\nP(F_0 = %d, F_1, F_2)\n", f_0);

        // top line
        printf("\n----------------");
        for (uint8_t f = 0; f < (this->f_max + 1); f++)
            printf("-------------");

        // line with values of f
        printf("\n[(F_1, F_2)] : |");
        for (uint8_t f = 0; f < (this->f_max + 1); f++)
            printf(" %-11d|", f);

        printf("\n----------------");
        for (uint8_t f = 0; f < (this->f_max + 1); f++)
            printf("-------------");

        for (uint8_t f_1 = 0; f_1 < (this->f_max + 1); f_1++) {

            iface_pivots[1] = f_1;

            printf("\n %-14d|", f_1);

            for (uint8_t f_2 = 0; f_2 < (this->f_max + 1); f_2++) {

                iface_pivots[2] = f_2;  

                printf(" %-.5LE|", this->get_joint_f_pmf(iface_pivots));

                total += this->get_joint_f_pmf(iface_pivots);
            }
        }

        printf("\n----------------");
        for (uint8_t f = 0; f < (this->f_max + 1); f++)
            printf("-------------");

        printf("\n");
    }

    printf("\nSUM(P(F_0, F_1, F_2)) = %-.5LE\n", total);
}

__float080 * RID_Router::init_lpm_iface_pmf() {

    // initialize 1 x (iface_num + 2) matrix for I's pmf
    this->lpm_iface_pmf = (__float080 *) calloc(this->iface_num + 2, sizeof(__float080));

    return this->lpm_iface_pmf;
}

void RID_Router::set_lpm_iface_pmf(int iface, __float080 value) {
    this->lpm_iface_pmf[iface + 2] = value;
}

__float080 RID_Router::get_lpm_iface_pmf(int iface) {
    return this->lpm_iface_pmf[iface + 2];
}

void RID_Router::print_lpm_iface_pmf() {

    printf("\n----------------");
    for (int i = IFACE_MULTI_HITS; i < (this->iface_num + 1); i++)
        printf("-------------");

    printf("\n[I]          : |");

    for (int i = -2; i < this->iface_num; i++)
        printf(" %-11d|", i);

    printf(" SUM P(I)   |");

    printf("\n----------------");
    for (int i = IFACE_MULTI_HITS; i < (this->iface_num + 1); i++)
        printf("-------------");

    __float080 _prob = 0.0;
    __float080 _prob_total = 0.0;

    printf("\n[P(I = i)] :   |");

    for (int i = IFACE_MULTI_HITS; i < this->iface_num; i++) {

        _prob = this->get_lpm_iface_pmf(i);
        _prob_total += _prob;

        printf(" %-.5LE|", _prob);
    }

    printf(" %-.5LE|", _prob_total);

    printf("\n----------------");
    for (int i = IFACE_MULTI_HITS; i < (this->iface_num + 1); i++)
        printf("-------------");

    printf("\n");
}

int RID_Router::calc_f_pmf(
    uint8_t request_size, 
    uint8_t ingress_iface,
    int * tp_sizes, 
    __float080 * f_r_distribution) {

    // placeholder for P(F_i = f)
    __float080 log_prob_Li_f;
    __float080 prob_Li_f;

    for (uint8_t i = IFACE_LOCAL; i < this->iface_num; i++) {

        // if iface i isn't initialized, abort
        if (this->fwd_table[i].iface == -1) {

            fprintf(stderr, "RID_Router::calc_f_pmf() : entry %d uninitialized. aborting.\n", i);
            return -1;
        }

        this->init_f_pmf(i);

        // let's assume that, according to LPM rules, we don't forward 
        // packets over the ingress iface (P(|F| = 0) = 1.0, and 0.0 for 
        // all f > 0)
        if (i == ingress_iface) {

            this->set_f_pmf(i, 0, 1.0);
            continue;
        }

        // get the 1 x f_max matrix with the probabilities of NOT 
        // having FPs among the entries of size |F| = f in iface i. 
        // these probabilities are given in log() form, since we'll 
        // be dealing with small numbers.
        __float080 * log_fp_rates = (__float080 *) calloc(this->f_max, sizeof(__float080));
        __float080 k = (log(2) * DEFAULT_M) / ((__float080) request_size);

        // printf("RID_Router::calc_f_pmf() : r[%d][%d]->fwd_table[%d].iface_proportion = %-.5LE\n", 
        //     this->height, this->width, i, this->fwd_table[i].iface_proportion);

        if (get_log_fp_rates(
            DEFAULT_M, (__float080) request_size, k, 
            this->fwd_table[i].iface_proportion, this->fwd_table[i].f_distribution, f_r_distribution, 
            log_fp_rates) < 0) {

            free(log_fp_rates);
            fprintf(stderr, "RID_Router::calc_f_pmf() : couldn't retrieve log FP rates\n");

            return -1;
        }

        // calculate P(F_i = f)
        for (int f = this->f_max; f >= 0; f--) {

            log_prob_Li_f = 0.0;
            prob_Li_f = 0.0;

            if (f < tp_sizes[i]) {

                // if f is smaller that the max. TP size, then according 
                // to LPM semantics, an entry with size f will never 
                // be chosen as the longest match: instead 
                // the max. TP entry is followed, not matter the 
                // iface it points to. thus, P(F_i = f) = 0.0
                this->set_f_pmf(i, f, 0.0);

            } else {

                // if f is larger than or equal to max. TP size, then 
                // that entry MAY be chosen as the longest match. f 
                // can be chosen as the longest match if the 
                // following events happen:
                // 1) no entry with |F| > f triggers a FP 
                // 2) at least 1 entry with |F| = f is chosen

                // note that event 2 can happen in 2 cases:
                // 2.1) f == max_tp_size and i == max_tp_iface : the tp 
                //      entry will be chosen for sure, P('event 2') = 1
                // 2.2) at least 1 entry with |F| = f triggers a FP

                // we assume events 1 and 2 are independent, 
                // therefore P(F_i = f) = P('event 1') * P('event 2')
                for (int _f = this->f_max; _f >= f; _f--) {

                    // event 1 : sum the logs of the probabilities of 
                    // not having FPs for |F| > f. this yields 
                    // P('event 1')
                    if (_f > f) {

                        log_prob_Li_f += log_fp_rates[_f - 1];

                    // event 2
                    } else {

                        // in the special case that f == max_tp_size 
                        // and i == max_tp_iface, then P('event 2') = 1
                        if (f == tp_sizes[i]) {

                            // P('event 1') * P('event 2')
                            prob_Li_f = exp(log_prob_Li_f) * (1.0);

                        } else {

                            prob_Li_f = exp(log_prob_Li_f) * (1.0 - exp(log_fp_rates[_f - 1]));
                        }
                    }
                }
            }

            // save P(F_i = f)
            this->set_f_pmf(i, f, prob_Li_f);
        }

        free(log_fp_rates);
    }

    return 0;
}

int RID_Router::calc_joint_f_pmf() {

    this->init_joint_f_pmf();

    // hold the pivots of the different iface arrays 
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int)); 

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;

    __float080 log_prob = 0.0;
    __float080 total_prob = 0.0;

    while (curr_i > -1) {

        // the inner loop is where all the work gets done
        if (curr_i == (this->iface_num - 1)) {

            for (int i = 0; i < (this->f_max + 1); i++) {

                // update the |F| size pivot for iface i 
                iface_pivots[this->iface_num - 1] = i;

                // calculate the log probability for the joint event 
                // encoded in iface_pivots
                log_prob = calc_joint_f_log_prob(iface_pivots);

                // set this value on the joint_f_pmf matrix kept by the 
                // path state variable 
                this->set_joint_f_pmf(iface_pivots, exp(log_prob));
                total_prob += exp(log_prob);
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

__float080 RID_Router::calc_cumulative_joint_f(uint8_t iface, uint8_t f) {

    __float080 cumulative_prob = 0.0;

    // hold the pivots of the different iface arrays 
    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int)); 

    int curr_i = (this->iface_num - 1);
    int curr_f = 0;

    // FIXME: this is sort of hacky, but, well... it works ?
    uint8_t going_up = 0x00;

    // this should never change for the rest of the computations
    iface_pivots[iface] = f;

    while (curr_i > -1) {

        // the inner loop is where all the work gets done
        if (curr_i == (this->iface_num - 1)) {

            if (curr_i != iface) {

                for (int i = 0; i < f; i++) {

                    // update the |F| size pivot for iface i 
                    iface_pivots[this->iface_num - 1] = i;
                    cumulative_prob += this->get_joint_f_pmf(iface_pivots);
                }

            } else {

                cumulative_prob += this->get_joint_f_pmf(iface_pivots);
                going_up = 0x00;
            }
        }

        // extract the current pivot value for the curr_i level
        if (curr_i == iface && going_up) {
            
            curr_f = 0;

        } else {

            curr_f = iface_pivots[curr_i];
        }

        if (++curr_f < f) {

            // we're at some level curr_i, and we need to increment 
            // it's pivot counter (to curr_f)
            if (curr_i != iface)
                iface_pivots[curr_i] = curr_f;

            // after updating iface_pivots[curr_i], go back down 
            // to the level below.
            curr_i++;
            going_up = 0xFF;

        } else {

            // we've completed the iterations for the curr_i level. we 
            // now: 
            //  * reset iface_pivots[curr_i] to 0
            if (curr_i != iface) {
                curr_f = -1;
                iface_pivots[curr_i] = curr_f;
            }

            //  * go up one level (i.e. to curr_i - 1) to increment 
            //      iface_pivots[curr_i - 1]
            curr_i--;
            going_up = 0x00;
        }
    }

    free(iface_pivots);

    return cumulative_prob;
}

__float080 RID_Router::calc_no_match_prob() {

    int * iface_pivots = (int *) calloc(this->iface_num, sizeof(int));

    __float080 drop_probability = this->get_joint_f_pmf(iface_pivots);

    free(iface_pivots);

    return drop_probability;
}

int RID_Router::calc_lpm_iface_pmf() {

    __float080 _prob_exclusive;
    __float080 _prob_forward = 0.0;
    __float080 _prob_aux = 0.0;

    this->init_lpm_iface_pmf();

    for (uint8_t i = 0; i < this->iface_num; i++) {

        _prob_exclusive = 0.0;

        for (uint8_t f = f_max; f > 0; f--) {

            _prob_exclusive += calc_cumulative_joint_f(i, f);
        }

        _prob_forward += _prob_exclusive;
        this->set_lpm_iface_pmf(i, _prob_exclusive);
    }

    this->set_lpm_iface_pmf(IFACE_NO_HITS, this->calc_no_match_prob());

    // FIXME: THIS IS SUPER HACK-ISH AND YOU SHOULD FEEL BAD...

    // i had to resort to this because we were getting negative probability 
    // values, with very small exponents (e.g. -17), which *CAN* be explained 
    // by approx. errors. that's the basis for the hack below. however, i 
    // haven't confirmed this yet...
    _prob_aux = (1.0 - (this->calc_no_match_prob() + _prob_forward));

    if (fabs(_prob_aux) < pow(10.0, -15))
        _prob_aux = 0.0;

    this->set_lpm_iface_pmf(IFACE_MULTI_HITS, _prob_aux);

    return 0;
}

__float080 RID_Router::calc_joint_f_log_prob(int * iface_pivots) {

    __float080 log_prob = 0.0;

    for (uint8_t i = 0; i < this->iface_num; i++)
        log_prob += log(this->get_f_pmf(i, iface_pivots[i])); 

    return log_prob;
}
