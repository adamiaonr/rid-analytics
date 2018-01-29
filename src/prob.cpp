#include "prob.h"

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

void Prob::calc_log_prob_fp_neq(
    fp_data iface_fp_data,
    std::vector<__float080> & log_prob_fp_neq) {

    // as a precaution, reset log_prob_fp_neq to 0.0
    std::fill(log_prob_fp_neq.begin(), log_prob_fp_neq.end(), 0.0);

    for (uint8_t f = 0; f < log_prob_fp_neq.size(); f++) {
        // total nr. of entries w/ size f
        __float080 num_entries = iface_fp_data.num_entries * iface_fp_data.f_distr[f];
        // std::cout << "Prob::calc_log_prob_fp_neq() : num_entries[" << f << "] = " 
        //     << num_entries << " (out of " << iface_fp_data.num_entries << ")" << std::endl;

        if (f == 0) {
            
            log_prob_fp_neq[f] += num_entries * log(1.0 - fp_rate(this->m, this->n, this->k, (__float080) (f + 1)));

        } else {

            // for entries w/ size f, we can have entries w/ |f\r| in the 
            // interval [1, f]. remember: the fp rates vary w/ |f\r|.
            __float080 subtotal_f_r = 0.0;
            // % of entries (of those w/ size = f) which have |f\r| <= f
            for (uint8_t j = 0; j <= f; j++) 
                subtotal_f_r += iface_fp_data.f_r_distr[j];

            for (uint8_t j = 0; j <= f; j++) {
                // total nr. of entries w/ size f and |f\r| = j
                __float080 num_fr_entries = (__float080) num_entries * (iface_fp_data.f_r_distr[j] / subtotal_f_r);
                // fp rate for size f, |f\r| = j
                log_prob_fp_neq[f] += num_fr_entries * log(1.0 - fp_rate(this->m, this->n, this->k, (__float080) (j + 1)));
            }
        }
    }
}

int Prob::calc_lm_prob(fp_data iface_fp_data, uint8_t i, bool anti) {

    // decide if we're calculating the iface or 'anti-iface' largest match 
    // probabilities
    std::vector<std::vector<std::vector<__float080> > > * lmp;
    if (!anti) lmp = &(this->lm_prob);
    else lmp = &(this->anti_lm_prob);

    std::cout << "Prob::calc_lm_prob(iface) : [INFO]: " 
        << "\t\n [n, m, k] : " << this->n << ", " << this->m << ", " << this->k
        << "\t\n [iface] : " << ((int) i)
        << "\t\n [# entries] : " << iface_fp_data.num_entries
        << "\t\n [lmp size] : " << (*lmp).size()
        << "\t\n [f_distr] : ";

    for(uint8_t f = 0; f < this->n; f++) 
        std::cout << "[" << (int) f << "] = " << iface_fp_data.f_distr[f] << ", ";
    std::cout << std::endl;

    // for (uint8_t j = 0; j < (*lmp)[i].size(); j++) {
    //     for (uint8_t k = 0; k < (*lmp)[i][j].size(); k++) {                  

    //         std::cout << "Prob::Prob() : [INFO] lm_prob[" << (int) i << "][" << (int) j << "][" 
    //             << (int) k << "] = " << (*lmp)[i][j][k] << std::endl;
    //     }
    // }

    // // if iface isn't initialized, abort
    // if (router->fwd_table[i].iface == -1) {
    //     std::cerr << "Prob::calc_lm_prob(iface) : [ERROR] iface " << i << " uninitialized. skipping." << std::endl;
    //     return -1;
    // }

    // aux. arrays used in fp rate calculation:
    //  - prob of not having an fp match of size f
    std::vector<__float080> log_prob_fp_neq(n, 0.0);
    //  - prob of fp match smaller than or equal to f
    std::vector<__float080> log_prob_fp_smeq(n + 1, 0.0);

    // if there are no entries for this iface, a match is impossible for 
    // this [i][ptree_size][f] combination
    // as such, match size 0 has prob 1.0 of happening
    if (iface_fp_data.num_entries < 1.0) {

        for (int fptree_size = (int) this->n; fptree_size >= 0; fptree_size--) {
            // std::cout << "Prob::calc_lm_prob(iface) : [INFO]: " 
            //     << "\t\n [fptree_size] : " << ((int) fptree_size) << std::endl;
            (*lmp)[i][fptree_size][0] = 1.0;
        }

        return 0;
    }

    // fill log_prob_fp_neq w/ prob of *only* having fp matches larger than f
    // FIXME : get a fp_data struct
    calc_log_prob_fp_neq(iface_fp_data, log_prob_fp_neq);

    // prob of fp match smaller than or equal to f 
    for (int f = (int) (this->n - 1); f >= 0; f--)
        log_prob_fp_smeq[f] = log_prob_fp_smeq[f + 1] + log_prob_fp_neq[f];

    // std::cout << "Prob::calc_lm_prob(iface) : [INFO] probabilities : \n\tp(fp > x)  : ";
    // for(uint8_t f = 0; f < this->n; f++) 
    //     std::cout << "[" << (int) (f) << "] = " << log_prob_fp_neq[f] << ", ";
    // std::cout << "\n\tp(fp <= x) : ";
    // for(uint8_t f = 0; f < (this->n + 1); f++) 
    //     std::cout << "[" << (int) f << "] = " << log_prob_fp_smeq[f] << ", ";
    // std::cout << std::endl;

    for (int fptree_size = (int) this->n; fptree_size >= 0; fptree_size--) {
        for (int f = (int) this->n; f >= 0; f--) {

            if (f < iface_fp_data.tp_size) {

                // if f is smaller than the max. true positive match 
                // then size f *WON'T BE* the largest match for sure
                (*lmp)[i][fptree_size][f] = 0.0;

            } else {

                // 2 sub-cases (deterministic):
                //  1) TP exists and its size is f
                //  2) no TP exists for f
                if (f == iface_fp_data.tp_size) {

                    // sub-case 1) TP exists and its size is f

                    // if iface is in the ptree and f < ptree_size : 
                    // since ptree_size < iface_fp_data.tp_size, f will never 
                    // be the largest match
                    if (f < fptree_size) {

                        // FIXME: this seems hack-ish : here i check if 
                        // iface is associated with any entries. why set 
                        // it to 1.0 if the iface has no entries? is this  
                        // even possible to reach?
                        (*lmp)[i][fptree_size][f] = 
                            (__float080) ((iface_fp_data.num_entries > 0) ? 0.0 : 1.0);

                    } else {

                        // if f >= ptree_size, the largest match will *AT LEAST* be f.
                        (*lmp)[i][fptree_size][f] = exp(log_prob_fp_smeq[f]);
                    }

                } else {

                    // sub-case 2) no TP exists for f
                    if (f > 0) {

                        // // if f is less than the prefix tree size, then 
                        // // any match will be larger than f
                        // if (f < ptree_size) {

                        //     this->(*lmp)[iface][ptree_size].(*lmp)[f] = 0.0;

                        // } else if (f == ptree_size) {
                        if (f == fptree_size) {

                            // the prob of this event can be calculated as:
                            //  * not having fp matches larger than f : log_prob_fp_smeq[f]
                            //  * AND having a fp match of size f : 1.0 - exp(log_prob_fp_neq[f - 1])
                            (*lmp)[i][fptree_size][f] = 
                                exp(log_prob_fp_smeq[f]) 
                                * (1.0 - exp(log_prob_fp_neq[f - 1])) 
                                * ((iface_fp_data.num_entries > 0) ? 1.0 : 0.0);

                        } else {

                            // if the iface is potentially on ptree, then it can't have 
                            // a match larger than the prefix tree. if 
                            // the iface is def. not on ptree, it can be the start of 
                            // a new prefix tree, and as such fp matches 
                            // larger than f are possible.
                            if (iface_fp_data.on_fptree && fptree_size > 0) {

                                (*lmp)[i][fptree_size][f] = 0.0;
                            
                            } else {

                                // the prob of this event can be calculated as:
                                //  * not having fps larger than f : log_prob_fp_smeq[f]
                                //  * AND having a fp of size f : 1.0 - exp(log_prob_fp_neq[f - 1])
                                (*lmp)[i][fptree_size][f] = 
                                    exp(log_prob_fp_smeq[f]) 
                                    * (1.0 - exp(log_prob_fp_neq[f - 1])) 
                                    * ((iface_fp_data.num_entries > 0) ? 1.0 : 0.0);
                            }
                        }

                    } else {

                        if (f < fptree_size) {

                            (*lmp)[i][fptree_size][f] = 
                                (__float080) ((iface_fp_data.num_entries > 0) ? 0.0 : 1.0);

                        } else {

                            (*lmp)[i][fptree_size][f] = 
                                exp(log_prob_fp_smeq[f]);
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int Prob::calc_lm_probs(std::vector<std::vector<fp_data> > * iface_fp_data) {

    if (this->has_lm_prob) {
        std::cout << "Prob::calc_lm_probs() : [INFO] pre-calculated lm_prob table available. skipping." << std::endl;
        return 0;
    }

    for (uint8_t i = 0; i < this->iface_num; i++) {

        calc_lm_prob((*iface_fp_data)[0][i], i);
        calc_lm_prob((*iface_fp_data)[1][i], i, true);

        print_lm_prob(i);
        print_lm_prob(i, true);
    }

    // mark lm probs as calculated
    this->has_lm_prob = true;

    return 0;
}

// objective: calc P('iface i on fp tree of size p'), for each iface i and fp 
//            tree size p
// 
// algorithm:
//  - for each [i][p] pair:
//      - [i][p = 0] : represents the case of iface i not on any *existing* fp tree,
//                     i.e. a fp tree that comes from a previous router
//          - in this case, we assume P([i][0]) = (1.0 / (num_valid_ifaces)), i.e. 
//            a uniform distribution of P(p = 0) per iface
//
//      - [i][p > 0] : iface i *may* be on a fp tree of size p.
//          - if iface_fp_data[i].on_fptree[p] == false, P([i][p]) = 0.0
//          - else, P([i][p]) = iface_fp_data[i].entry_prop * iface_fp_data[i].f_ditr[p] / (f_prop[p]), 
//            i.e. the proportion of entries of size p held by iface i
//
// complexity: O(N*I)
//
int Prob::calc_iface_on_fptree_probs(
    std::vector<std::vector<fp_data> > * iface_fp_data,
    std::vector<__float080> * in_fptree_prob) {

    // calc number of valid ifaces
    __float080 num_valid_ifaces = 0.0;
    for (uint8_t i = 0; i < this->iface_num; i++) {
        if (!(*iface_fp_data)[0][i].is_blocked) num_valid_ifaces++;
    }

    // calc the entry size proportions, i.e. share of entries of size f among 
    // all entries in the router.
    std::vector<__float080> f_prop(this->n, 0.0);
    for (uint8_t f = 0; f < this->n; f++) {
        for (uint8_t i = 0; i < this->iface_num; i++) {
            f_prop[f] += (*iface_fp_data)[0][i].entry_prop * (*iface_fp_data)[0][i].f_distr[f];
        }
    }

    // finally, calc the actual iface on fp tree probabilities
    for (uint8_t i = 0; i < this->iface_num; i++) {

        if ((*iface_fp_data)[0][i].is_blocked) continue;

        // for P([i][p = 0])
        // this->iface_on_fptree_prob[i][0] = (*in_fptree_prob)[0] * (1.0 / num_valid_ifaces);
        this->iface_on_fptree_prob[i][0] = (*in_fptree_prob)[0];
        // for P([i][p > 0])
        for (uint8_t p = 1; p < (this->n + 1); p++) {

            if (!(*iface_fp_data)[0][i].on_fptree) break;

            this->iface_on_fptree_prob[i][p] = (*in_fptree_prob)[p] * (*iface_fp_data)[0][i].entry_prop * (*iface_fp_data)[0][i].f_distr[p - 1];
            this->iface_on_fptree_prob[i][p] /= f_prop[p - 1];
        }
    }

    return 0;
}

__float080 Prob::calc_iface_prob(uint8_t i, bool strict) {

    __float080 iface_prob = 0.0;
    for (uint8_t fptree_size = 0; fptree_size < (this->n + 1); fptree_size++) {

        __float080 cond_prob = 0.0;
        for (int f = this->n; f >= 0; f--) {
            // in 'strict' mode, we only consider cases for which L_i > L_{~i}
            // else, we consider L_i >= L_{~i}
            int m = 0;
            (strict ? m = f : m = (f + 1));

            // calc the conditional prob P(L_i >= L_{~i} | 'i on p')
            //  - note that lm_prob holds P(L_i = f | 'i on p')
            //  - note that anti_lm_prob holds P(L_{~i} = k | 'i on p')
            
            //  - as such, we do P(L_i >= L_{~i} | 'i on p') = 
            //      = sum^(m)_k [ P(L_i = f, L_{~i} = k | 'i on p') ]
            //      = sum^(m)_k [ P(L_i = f | 'i on p') * P(L_{~i} = k | 'i on p') }

            //  - note that we assume 'i on p' == '~i on 0', i.e. 
            //    if iface i is on the fp tree p, then ~i is not on any fp tree  
            for (int k = 0; k < m; k++)
                cond_prob += (this->lm_prob[i][fptree_size][f] * this->anti_lm_prob[i][0][k]);
        }

        // by def of joint prob : P(L_i >= L_{~i} , 'i on p') = P(L_i >= L_{~i} | 'i on p') * P('i on p')
        iface_prob += cond_prob * this->iface_on_fptree_prob[i][fptree_size];
    }

    return iface_prob;
}

// objective: calc P(i) = P(L_i >= L_{~i}), for each i, i.e. the probability of 
//            choosing iface i on a forwarding decision
//
// algorithm:
//
//  - for each iface i, calc the conditional probs: 
//      - P(L_i >= L_{~i} | 'i is on fptree p'), for all p
//
//      - the complexity of each of the steps above is O(I*N^3). 
//          - N is the RID req. size (or, alternatively, the nr. of 
//            elements which can be encoded in a bf). 
//          - I is the nr. of ifaces in the router
//          
//      - why?
//          - we compare each of the N different values of L_{i,p}, with 
//            L_{~i,0}. that makes O(N^2). 
//          - if we do the above for all p (which is equal to N), that makes 
//            O(N^3).
//          - finally, doing the above for all I ifaces, that makes O(I*N^3)
//          - while this may seem like a crazy complexity, not it is better 
//            than the previous O(N^(I + 1)).
//
//  - finally, marginalize P(L_i >= L_{~i} | 'i in p') above over all p, to get 
//    P(i) (for all ifaces i). the final complexity becomes O(I*N^4)
//
// complexity: O(I*N^4), I being the nr. of ifaces, N the req. size
//
int Prob::calc_iface_probs(
    std::vector<std::vector<fp_data> > * iface_fp_data,
    std::vector<__float080> * in_fptree_prob,
    std::vector<__float080> & iface_probs,
    bool strict) {

    if (calc_lm_probs(iface_fp_data) < 0) 
        return -1;

    if (calc_iface_on_fptree_probs(iface_fp_data, in_fptree_prob) < 0)
        return -1;

    for (uint8_t i = 0; i < this->iface_num; i++)
        iface_probs[i] = calc_iface_prob(i, strict);

    return 0;
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
void Prob::print_lm_prob(uint8_t iface, bool anti) {

    // decide if we're calculating the iface or 'anti-iface' largest match 
    // probabilities
    std::vector<std::vector<std::vector<__float080> > > * lmp;
    if (!anti) lmp = &(this->lm_prob);
    else lmp = &(this->anti_lm_prob);

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->n + 2); i++)
        printf("-------------");

    printf("\n[|L|]         : |");

    for (uint8_t i = 0; i < (this->n + 1); i++)
        printf(" %-11d|", i);

    if (!anti) printf(" SUM P(L_{i,p}) |");
    else printf(" SUM P(L_{~i,p}) |");

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->n + 2); i++)
        printf("-------------");

    __float080 cumulative_prob = 0.0;
    __float080 prob = 0.0;

    for (uint8_t fptree_size = 0; fptree_size < (this->n + 1); fptree_size++) {

        if (!anti) printf("\n[P(L_{%d,%d}) ] : |", iface, fptree_size);
        else printf("\n[P(L_{~%d,%d})] : |", iface, fptree_size);

        for (uint8_t f = 0; f < (this->n + 1); f++) {

            prob = (*lmp)[iface][fptree_size][f];
            printf(" %-.5LE|", prob);

            cumulative_prob += prob;
        }

        printf("     %-.5LE|", cumulative_prob);
        cumulative_prob = 0.0;
    }

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->n + 2); i++)
        printf("-------------");

    printf("\n");
}