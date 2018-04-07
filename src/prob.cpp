#include "prob.h"

__float080 fp_rate(__float080 m, __float080 n, __float080 k, __float080 c) {
    return pow((1.0 - exp(-((n / m) * k))), k * c);
}

int clear_lsb(int & number) {
    number = (number & (number - 1));
    return number;
}

int count_set_bits(int number) {

    int set_bits = 0;
    int _num = number;
    
    while(_num > 0) {
        clear_lsb(_num);
        set_bits++;
    }

    return set_bits;
}

void Prob::calc_log_prob_fp_neq(
    Prob::fp_data * iface_fp_data,
    std::vector<__float080> & log_prob_fp_neq,
    std::vector<__float080> lt_ratios) {

    // as a precaution, reset log_prob_fp_neq to 0.0
    std::fill(log_prob_fp_neq.begin(), log_prob_fp_neq.end(), 0.0);

    for (uint8_t f = 0; f < log_prob_fp_neq.size(); f++) {

        // total nr. of entries w/ size f
        __float080 num_entries = iface_fp_data->num_entries * iface_fp_data->f_distr[f] * lt_ratios[f];
        // std::cout << "Prob::calc_log_prob_fp_neq() : num_entries[" << f << "] = " 
        //     << num_entries << " (out of " << iface_fp_data->num_entries << ")" << std::endl;

        if (f == 0) {
            
            log_prob_fp_neq[f] += num_entries * log(1.0 - fp_rate(this->m, this->n, this->k, (__float080) (f + 1)));

        } else {

            // for entries w/ size f, we can have entries w/ |f\r| in the 
            // interval [1, f]. remember: the fp rates vary w/ |f\r|.
            __float080 subtotal_f_r = 0.0;
            // % of entries (of those w/ size = f) which have |f\r| <= f
            for (uint8_t j = 0; j <= f; j++) 
                subtotal_f_r += iface_fp_data->f_r_distr[j];

            for (uint8_t j = 0; j <= f; j++) {
                // total nr. of entries w/ size f and |f\r| = j
                __float080 num_fr_entries = (__float080) num_entries * (iface_fp_data->f_r_distr[j] / subtotal_f_r);
                // fp rate for size f, |f\r| = j
                log_prob_fp_neq[f] += num_fr_entries * log(1.0 - fp_rate(this->m, this->n, this->k, (__float080) (j + 1)));
            }
        }
    }
}

int count_shared_trees(
    int t,
    Prob::fp_data * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks) {

    int count = 0;
    for (unsigned int k = 0; k < (*tree_bitmasks)[t].size(); k++)
        count += count_set_bits((int) (((*tree_bitmasks)[t][k]) & (iface_fp_data->tree_bitmasks[t][k])));

    return count;
}

void local_tree_ratios(
    Prob::fp_data * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,
    std::vector<__float080> & lt_ratios) {

    std::map<int, std::vector<uint8_t> >::iterator itr;
    for (itr = (*tree_bitmasks).begin(); itr != (*tree_bitmasks).end(); itr++) {

        // isolate the tree size
        int t = itr->first;
        // skip 'union' of tree bitmasks
        if (t == 0) continue;
        // count trees & shared trees which pass through iface
        int trees = 0, shared_trees = 0;
        for (unsigned int k = 0; k < (*tree_bitmasks)[t].size(); k++) {
            trees += count_set_bits((int) iface_fp_data->tree_bitmasks[t][k]);
            shared_trees += count_set_bits((int) (((*tree_bitmasks)[t][k]) & (iface_fp_data->tree_bitmasks[t][k])));
        }

        if (trees > 0)
            lt_ratios[t - 1] = (__float080) (trees - shared_trees) / (__float080) (trees);
        std::cout << "\t [# trees] : " << trees
            << "\n\t [# shared trees] : " << shared_trees
            << "\n\t [lt ratios[" << (t) << "]] : " << lt_ratios[t - 1] << std::endl;
    }
}

void _calc_lm_prob(
    int f, int t,
    Prob::fp_data * iface_fp_data,
    std::vector<__float080> * log_prob_fp_neq,
    std::vector<__float080> * log_prob_fp_smeq,
    std::vector<__float080> * lmp) {

    // if f is smaller than the tp size or fp tree size t, then f *WON'T 
    // BE* the largest match for sure, so P(L_i = f) = 0.0
    if ((f < iface_fp_data->tp_size) || f < t) {

        (*lmp)[f] = 0.0;

    } else {

        // 2 sub-cases (deterministic):
        //  1) f == tp size
        //  2) f > tp size
        if (f == iface_fp_data->tp_size) {

            // 1) f == tp size
            // the prob of this event can be calculated as:
            //  * not having fp matches larger than f : log_prob_fp_smeq[f]
            (*lmp)[f] = exp((*log_prob_fp_smeq)[f]);

        } else {

            // 2) f > tp size
            if (f == t) {

                // the prob of this event can be calculated as:
                //  * not having fp matches larger than f : log_prob_fp_smeq[f]
                (*lmp)[f] = exp((*log_prob_fp_smeq)[f]);

            } else {

                if (iface_fp_data->on_fptree[t] && t > 0) {

                    (*lmp)[f] = 0.0;
                
                } else {

                    // the prob of this event can be calculated as:
                    //  * not having fps larger than f : log_prob_fp_smeq[f]
                    //  * AND having a fp of size f : 1.0 - exp(log_prob_fp_neq[f - 1])
                    (*lmp)[f] = 
                        exp((*log_prob_fp_smeq)[f]) 
                        * (1.0 - exp((*log_prob_fp_neq)[f - 1])) 
                        * ((iface_fp_data->num_entries > 0) ? 1.0 : 0.0);
                }
            }
        }
    }
}

int Prob::calc_lm_prob(
    uint8_t i, 
    uint8_t t,
    Prob::fp_data * iface_fp_data,
    std::vector<__float080> * lmp) {

    // aux. arrays used in fp rate calculation:
    //  - P(|FP| =/= f) = (1.0 - P(|FP| = f)) : prob of not having an fp match of size f
    std::vector<__float080> log_prob_fp_neq(n, 0.0);
    //  - P(|FP| <= f) : prob of fp match smaller than or equal to f
    std::vector<__float080> log_prob_fp_smeq(n + 1, 0.0);

    if (iface_fp_data->num_entries < 1.0) { (*lmp)[0] = 1.0; return 0; }

    // fill P(|FP| =/= f) = (1.0 - P(|FP| = f)), for all f
    std::vector<__float080> lt_ratios = std::vector<__float080>(n, 1.0);
    calc_log_prob_fp_neq(iface_fp_data, log_prob_fp_neq, lt_ratios);
    // calc P(|FP| <= f), for all f
    for (int f = (int) (this->n - 1); f >= 0; f--)
        log_prob_fp_smeq[f] = log_prob_fp_smeq[f + 1] + log_prob_fp_neq[f];

    for (int f = (int) this->n; f >= 0; f--) {

        _calc_lm_prob(
            f, t, 
            iface_fp_data, 
            &log_prob_fp_neq, &log_prob_fp_smeq,
            lmp);
    }

    return 0;
}

int Prob::calc_lm_prob(
    uint8_t i, 
    Prob::fp_data * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,
    bool iface_complement) {

    // either 'i' or '~i (complement)' largest match probs
    std::vector<std::vector<std::vector<__float080> > > * lmp;
    if (!iface_complement) lmp = &(this->lm_cond_prob);
    else lmp = &(this->lm_complement_cond_prob);

    // aux. arrays used in fp rate calculation:
    //  - P(|FP| =/= f) = (1.0 - P(|FP| = f)) : prob of not having an fp match of size f
    std::vector<__float080> log_prob_fp_neq(n, 0.0);
    //  - P(|FP| <= f) : prob of fp match smaller than or equal to f
    std::vector<__float080> log_prob_fp_smeq(n + 1, 0.0);

    // if there are no entries for this iface, a match is impossible
    // as such, match size f = 0 has prob 1.0 of happening
    if (iface_fp_data->num_entries < 1.0) {
        for (int t = (int) this->n; t >= 0; t--) (*lmp)[i][t][0] = 1.0;
        return 0;
    }

    // calc the ratio of trees which can origin FP trees locally at the router
    std::cout << "Prob::calc_lm_prob(iface) : [INFO] lt ratio iface  " << (iface_complement ? "~" : "") << (int) i << " :" << std::endl;
    // std::vector<__float080> lt_ratios(n, 0.0);
    // local_tree_ratios(iface_fp_data, tree_bitmasks, lt_ratios);
    std::vector<__float080> lt_ratios = std::vector<__float080>(n, 1.0);

    std::cout << "Prob::calc_lm_prob(iface) : [INFO] iface " << (iface_complement ? "~" : "") << ((int) i) << " stats :"
        << "\n\t [n, m, k] : " << this->n << ", " << this->m << ", " << this->k << ")";
        // << "\n\t [# entries] : " << iface_fp_data->num_entries << " (" << (iface_fp_data->num_entries * lt_ratio) << ")";
    //     << "\t\n [lmp size] : " << (*lmp).size()
    //     << "\t\n [f_distr] : ";
    // for(uint8_t f = 0; f < this->n; f++) 
    //     std::cout << "[" << (int) f << "] = " << iface_fp_data->f_distr[f] << ", ";
    std::cout << std::endl;

    // fill P(|FP| =/= f) = (1.0 - P(|FP| = f)), for all f
    calc_log_prob_fp_neq(iface_fp_data, log_prob_fp_neq, lt_ratios);
    // calc P(|FP| <= f), for all f
    for (int f = (int) (this->n - 1); f >= 0; f--)
        log_prob_fp_smeq[f] = log_prob_fp_smeq[f + 1] + log_prob_fp_neq[f];

    for (int t = (int) this->n; t >= 0; t--) {
        for (int f = (int) this->n; f >= 0; f--) {
            _calc_lm_prob(
                f, t, 
                iface_fp_data, 
                &log_prob_fp_neq, &log_prob_fp_smeq,
                &((*lmp)[i][t]));
        }
    }

    return 0;
}

int Prob::calc_lm_probs(
    std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks) {

    for (uint8_t i = 0; i < this->iface_num; i++) {

        // time complexity : O(N^2 * I)
        calc_lm_prob(i, (*iface_fp_data)[0][i], tree_bitmasks);       // iface i
        // time complexity : O(N^2 * I)
        calc_lm_prob(i, (*iface_fp_data)[1][i], tree_bitmasks, true); // iface ~i

        print_lm_prob(i);
        print_lm_prob(i, true);
    }

    return 0;
}

// objective : calc P(T_{out,i} = t), for each iface i and all fp tree sizes t
//
// the objective is calculating the probability of the event 'the largest fp tree 
// size a request going out over iface i belongs to equals t', or P(T_{out,i} = t) in 
// short.
// 
// algorithm : 
//  - P(T_{out,i} = 0):
//      - iff iface_fp_data[i].tp_size > 0
//      - P(T_{out,i} = 0) = P(|FP| = 0)
//
//  - P(T_{out,i} > t) is influenced by 2 events:
//      1) a fp tree of size t - to which the request is stuck - continues on 
//         iface i
//      2) the request falls on a `fresh' fp tree of size t. 
//         this can happen due to 2 sub-events:
//          - 2.1 : request got into the router w/ no association w/ fp tree (P(T_in = 0))
//          - 2.2 : request got into the router assoc. w/ fp tree, but iface i 
//                  is not assoc. w/ fp tree
//
// time complexity : 
//  - O(N) : 
//      - for each iface i (this function)
//  - O(N * I) : 
//      - for all ifaces
int Prob::calc_out_fptree_prob(
    uint8_t i,
    Prob::fp_data * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,
    std::vector<__float080> * in_fptree_prob,
    std::vector<std::vector<__float080> > & out_fptree_probs) {

    // aux. arrays used in fp rate calculation:
    //  - P(|FP| =/= f) = (1.0 - P(|FP| = f)) : prob of not having an fp match of size f
    std::vector<__float080> log_prob_fp_neq(n, 0.0);
    //  - P(|FP| <= f) : prob of fp match smaller than or equal to f
    std::vector<__float080> log_prob_fp_smeq(n + 1, 0.0);

    // init all probs to 0.0
    for (uint8_t t = 0; t <= this->n; t++)
        out_fptree_probs[i][t] = 0.0;

    // we add the influence of events 2.1 and 2.2
    //  - k = 0 : event 2.1
    //  - k = 1 : event 2.2
    std::vector<__float080> lt_ratios;

    for (int k = 0; k < 2; k++) {

        __float080 subevent_prob = 0.0;

        if (k == 0) {
        
            // 2.1) req. got into router w/ no assoc. w/ fp tree
            subevent_prob = (*in_fptree_prob)[0];
            // all entries count for fp prob calculation
            lt_ratios = std::vector<__float080>(n, 1.0);

        } else {

            // 2.2) req. assoc. fp tree, but not iface i
            subevent_prob = (1.0 - (*in_fptree_prob)[0]) * (this->fptree_prob[0][i][0]);
            // only account w/ % of entries coming from local trees
            lt_ratios = std::vector<__float080>(n, 1.0);
            // local_tree_ratios(iface_fp_data, tree_bitmasks, lt_ratios);
        }

        calc_log_prob_fp_neq(iface_fp_data, log_prob_fp_neq, lt_ratios);
        for (int f = (int) (this->n - 1); f >= 0; f--)
            log_prob_fp_smeq[f] = log_prob_fp_smeq[f + 1] + log_prob_fp_neq[f];

        // P(T_{out,i} = 0)
        if (iface_fp_data->tp_size > 0) {

            // the req. can advance to the next router w/ t = 0 w/ prob P(|FP| = 0)
            // this can happen if:
            //  - iface i is not assoc. to fp tree
            //      - cases 2.1 or 2.2 above
            //  - P(FP = 0)
            out_fptree_probs[i][0] += 
                subevent_prob                   // i not associated w/ fp tree
                * exp(log_prob_fp_smeq[0]);     // not have a FP > 0 (i.e. having only FP = 0)

        } else {

            // there's no way of having a req. advance to the next router
            // without at least 1 fp match
            out_fptree_probs[i][0] = 0.0;
        }

        // P(T_{out,i} > 0)
        for (uint8_t t = 1; t <= this->n; t++) {

            out_fptree_probs[i][t] +=
                subevent_prob
                * exp(log_prob_fp_smeq[t])              // not have a FP > t
                * (1.0 - exp(log_prob_fp_neq[t - 1]));  // have a FP of size t
        }
    }

    // contribution of event 1
    for (uint8_t t = 1; t <= this->n; t++)
        out_fptree_probs[i][t] += (1.0 - (*in_fptree_prob)[0]) * this->fptree_prob[0][i][t];

    return 0;
}

// objective : calc P(T_i = t), for each iface i and all fp tree sizes t
// 
// the objective is calculating the probability of the events 
// 'largest fp tree size iface i belongs to is t', P(T_i = t) for short.
// this corresponds to transforming P(T_{in} = t) - or in_fptree_prob[] - 
// in P(T_i = t).
//
// this influences the calculation of P(L_i = l), P(L_i > L_{~i}) and 
// as ultimately P(iface = i).
// this is the most ambiguous part of the probability framework, and possibly 
// the part which is 'wrong'.
//
// algorithm : we have several options:
//  - weighted ifaces method:
//      - P(T_i = t) = w[i][t] * P(T_{in} = t)
//      - SUM [all i][all t] { P(T_i = t) } = SUM [all t] { P(T_{in} = t) }
//      - w[i][t] = e[i][t] * tr[i]
//          - e[i][t] : (# entries size t on iface i) / (# entries in router)
//          - tr[i]   : (# of trees in iface i shared w/ prev routers) / (# of trees shared w/ prev routers)
//      - rationale:
//          - <please fill this>
//
//  - other methods ?
//
// time complexity : 
//  - O(N * I)
int Prob::calc_fptree_probs(
    std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,
    std::vector<__float080> * in_fptree_prob) {

    // (1) calculate str[i][t] : shared tree ratios for iface i, size t
    __float080 str_total = 0.0;
    std::vector< std::vector<__float080> > str(this->iface_num, std::vector<__float080>(this->n + 1, 0.0));

    for (uint8_t i = 0; i < this->iface_num; i++) {
        // don't consider blocked ifaces
        if ((*iface_fp_data)[0][i]->is_blocked) continue;

        std::map<int, std::vector<uint8_t> >::iterator itr;
        for (itr = (*tree_bitmasks).begin(); itr != (*tree_bitmasks).end(); itr++) {

            int t = itr->first;
            if (t == 0) continue;

            // count shared trees w/ bitmask on iface i
            str[i][t] = (__float080) count_shared_trees(t, (*iface_fp_data)[0][i], tree_bitmasks);
            // keep track of total # of shared trees
            str_total += str[i][t];
        }
    }

    std::cout << "Prob::calc_fptree_probs() : shared tree ratios :" << std::endl;
    for (uint8_t i = 0; i < this->iface_num; i++)
        for (uint8_t t = 1; t < (this->n + 1); t++)
            std::cout << "\tstr[" << (int) i << "][" << (int) t << "] = " << str[i][t] / ((str_total > 0.0) ? str_total : 1.0) << std::endl;

    // (2) calc P(T_i = t), for all {i,t} pairs
    for (uint8_t i = 0; i < this->iface_num; i++) {
        
        // if iface i is blocked, don't bother going on...
        if ((*iface_fp_data)[0][i]->is_blocked) continue;
        // uint8_t max_tree = 0;
        __float080 cumulative_fptree_prob[2] = {0.0, 0.0};
        for (uint8_t t = 1; t < (this->n + 1); t++) {

            // P(T_i > 0) = srt[i][t] * P(T_in = t)
            this->fptree_prob[0][i][t] = (str[i][t] / ((str_total > 0.0) ? str_total : 1.0)) * (*in_fptree_prob)[t];
            // P(T_{~i} > 0) = P(T_in = t) - (srt[i][t] * P(T_in = t))
            this->fptree_prob[1][i][t] = ((*in_fptree_prob)[t] - this->fptree_prob[0][i][t]);

            // SUM {for all t > 0} [ srt[i][t] * P(T_in = t) ]
            //  - will be used in calc. of P(T_i = 0, B)
            cumulative_fptree_prob[0] += this->fptree_prob[0][i][t];
            cumulative_fptree_prob[1] += this->fptree_prob[1][i][t];
        }

        // P(T_i = 0) = 1.0 - SUM {for all t > 0} [ srt[i][t] * P(T_in = t) ]
        this->fptree_prob[0][i][0] = 1.0 - cumulative_fptree_prob[0];
        this->fptree_prob[1][i][0] = (__float080) ((float) 1.0 - (float) cumulative_fptree_prob[1]);

        std::cout << "Prob::calc_fptree_probs() : P(T[" << (int) i << "] = 0) = 1.0 - " 
            << cumulative_fptree_prob[0] << " = " << (this->fptree_prob[0][i][0]) << std::endl;
        std::cout << "Prob::calc_fptree_probs() : P(T[~" << (int) i << "] = 0) = 1.0 - " 
            << cumulative_fptree_prob[1] << " = " << (this->fptree_prob[1][i][0]) << std::endl;
    }

    std::cout << "Prob::calc_fptree_probs() : P(T_i = t) :" << std::endl;
    for (uint8_t i = 0; i < this->iface_num; i++)
        for (uint8_t t = 0; t < (this->n + 1); t++)
            std::cout << "\tP(T[" << (int) i << "] = " << (int) t << "]) = " << this->fptree_prob[0][i][t] << std::endl;

    std::cout << "Prob::calc_fptree_probs() : P(T_{~i} = t) :" << std::endl;
    for (uint8_t i = 0; i < this->iface_num; i++)
        for (uint8_t t = 0; t < (this->n + 1); t++)
            std::cout << "\tP(T[~" << (int) i << "] = " << (int) t << "]) = " << this->fptree_prob[1][i][t] << std::endl;

    return 0;
}

// objective: calc P(iface = i) for an iface i, i.e. the probability of 
//            choosing iface i on a forwarding decision
//
// algorithm: 
//  - to get P(iface i), we must go through a series of probability math... 
//    we explain it briefly below.
//
//  - for each iface i, we define that P(iface = i) = P(L_i >= L_{~i})
//      - why? : if the longest match on iface i is larger than (or equal) 
//               to the match sizes on its complement (i.e. all ifaces
//               other than i), then iface i will be chosen for 
//               forwarding.
//
//      - P(L_i >= L_{~i}) = ?
//          - P(L_i = l) = SUM[ALL P](P(L_i = l, T_i = t)) = SUM[ALL P]( P(L_i = l | T_i = t) * P(T_i = t) )
//              - we can do this because the events T_i = t: 
//                  - are mutually disjoint
//                  - their union is the entire sample space P_i
//          - same as above - with appropriate replacement of i by '~i' - for P(L_{~i} = l)
//          - as such we can calculate P(L_i >= L_{~i}) as:
//              - P(L_i >= L_{~i}) = SUM [l in [1:N]] { SUM [k in [0:l]] { P(L_i = l) * P(L_{~i} = k) } }
//
//      - how about time complexities?
//          - P(L_i = l) is O(N^2)
//          - P(L_i >= L_{~i}) is O(N^2)
//          - as such, for all ifaces i, complexity if O(N^2 * I)
//          - while this may seem like a crazy complexity, it's better 
//            than the previous O(N^(I)) (of course, this won't go into the 
//            paper, but just putting it out there...).
//
// time complexity: 
//  - O(N^2) : for any iface i
//  - O((N^2) * I) : for all ifaces
//
int Prob::calc_iface_prob(
    uint8_t i, 
    std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
    std::vector<__float080> * in_fptree_prob,
    std::vector<__float080> & iface_prob) {

    //  - P(L_i = l) and P(L_{~i} = l)
    std::vector<std::vector<__float080> > lm_marg_prob((this->iface_num), std::vector<__float080> ((this->n + 1), 0.0));
    std::vector<std::vector<__float080> > lm_complement_marg_prob((this->iface_num), std::vector<__float080> ((this->n + 1), 0.0));
    //  - P(L_i | T_i = 0, A) and P(L_~i | T_~i = 0, A)
    std::vector<std::vector<__float080> > event_a(2, std::vector<__float080> ((this->n + 1), 0.0));
    calc_lm_prob(i, 0, (*iface_fp_data)[0][i], &event_a[0]);
    calc_lm_prob(i, 0, (*iface_fp_data)[1][i], &event_a[1]);

    // 1) calc P(L_i = l) and P(L_{~i} = l), for all l
    for (int f = this->n; f >= 0; f--) {

        // P(L_i, T_i = 0) is a special case, influenced by 2 mutually exclusive 
        // events, say A and B:
        //
        //  - A : a request comes into the router while *not* associated w/ 
        //        a tree of any size t. 
        //        in this case, none of the ifaces can be associated w/ 
        //        a fp tree of any size, yielding:
        //          - P(T_i = 0) = 1.0
        //          - P(T_i > 0) = 0.0
        //
        //        this can happen iff a TP match and no FP matches happened 
        //        in the previous router.
        //        
        //  - B : say the req. comes into the router associated w/ a tree of some size t (t > 0).
        //        in this this case iface i can either be associated or not 
        //        associated w/ a fp tree of size t > 0, yielding:
        //          - P(T_i = 0) = 1.0 - SUM {for all t > 0} [ srt[i][t] * P(T_in = t) ]
        //          - P(T_i > 0) = srt[i][t] * P(T_in = t)
        //         
        // finally, P(L_i, T_i = 0) = 
        //      P(T_in = 0) * (P(L_i | T_i = 0, A) * P(T_i = 0, A))         // event A
        //      + (1 - P(T_in = 0)) * (P(L_i | T_i = 0, B) * T_i = 0, B)    // event B

        // - P(L_i, T_i = 0)
        lm_marg_prob[i][f] += 
            ((*in_fptree_prob)[0] * event_a[0][f] * ((*iface_fp_data)[0][i]->is_blocked ? 0.0 : 1.0))
            + ((1.0 - (*in_fptree_prob)[0]) * (this->lm_cond_prob[i][0][f] * this->fptree_prob[0][i][0]));
        // - P(L_~i, T_~i = 0)
        lm_complement_marg_prob[i][f] += 
            ((*in_fptree_prob)[0] * event_a[1][f])
            + ((1.0 - (*in_fptree_prob)[0]) * (this->lm_complement_cond_prob[i][0][f] * this->fptree_prob[1][i][0]));

        for (uint8_t t = 1; t < (this->n + 1); t++) {

            // - P(L_i = l, T_i = t)
            __float080 joint_prob = (1.0 - (*in_fptree_prob)[0]) * (this->lm_cond_prob[i][t][f] * this->fptree_prob[0][i][t]);
            // add joint P(L_i = l, T_i = t) to marginal P(L_i = l)
            lm_marg_prob[i][f] += joint_prob;

            // print if P(L_i = l, T_i = t) > 0.0
            if (joint_prob > 0.0) {
                std::cout << "\t[" << (int) t << "] +" 
                    << (joint_prob) << " (" << this->lm_cond_prob[i][t][f] << " x " << this->fptree_prob[0][i][t] << ")"
                    << " to P(L_" << (int) i << " = " << (int) f << ") = " << lm_marg_prob[i][f] << std::endl;
            }

            // - P(L_{~i} = l, T_{~i} = t)
            joint_prob = 
                    (1.0 - (*in_fptree_prob)[0])
                    * this->lm_complement_cond_prob[i][t][f]    // P(L_{~i} = l | T_{~i} = t)
                    * this->fptree_prob[1][i][t];               // P(T_{~i} = t)
            // add the joint P(L_{~i} = l, T_{~i} = t) to marginal P(L_{~i} = l)
            lm_complement_marg_prob[i][f] += joint_prob;
            // print if P(L_{~i} = l, T_{~i} = t) > 0.0
            if (joint_prob > 0.0) {
                std::cout << "\t[" << (int) t << "] +" 
                    << (joint_prob) << " (" << this->lm_complement_cond_prob[i][t][f] << " x " << this->fptree_prob[1][i][t] << ")"
                    << " to P(L_{~" << (int) i << "} = " << f << ") = " << lm_complement_marg_prob[i][f] << std::endl;
            }
        }
    }

    std::cout << "Prob::calc_iface_prob() : P(L_" << (int) i << " = l) :" << std::endl; 
    for (int f = 0; f < this->n + 1; f++)
        std::cout << "\tP(L_" << (int) i << " = " << (int) f << ") = " << lm_marg_prob[i][f] << std::endl;

    std::cout << "Prob::calc_iface_prob() : P(L_{~" << (int) i << "} = l) :" << std::endl; 
    for (int f = 0; f < this->n + 1; f++)
        std::cout << "\tP(L_{~" << (int) i << "} = " << (int) f << ") = " << lm_complement_marg_prob[i][f] << std::endl;


    // 2) calc P(L_i >= L_{~i})
    for (int f = this->n; f >= 0; f--) {

        __float080 prob_strict = 0.0;
        __float080 prob_total = 0.0;
        for (uint8_t k = 0; k < f; k++)
            prob_strict += (lm_marg_prob[i][f] * lm_complement_marg_prob[i][k]);

        prob_total = prob_strict + (lm_marg_prob[i][f] * lm_complement_marg_prob[i][f]);

        iface_prob[0] += prob_strict;
        iface_prob[1] += prob_total;
    }

    return 0;
}

// objective : calc all probs necessary for rid analytics:
//              - P(event = e)
//              - P(iface = i)
//              - P(longest fptree size = t)
//
// algorithm :
//  - check each of the functions for the respective algos
//
// time complexity :
//  - overall, O()
int Prob::calc_probs(
    std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
    std::map<int, std::vector<uint8_t> > * tree_bitmasks,    // tree bitmasks from prev router
    __float080 in_prob,
    std::vector<__float080> * in_fptree_prob,
    std::vector<std::vector<__float080> > & iface_probs,
    std::vector<std::vector<__float080> > & out_fptree_probs) {

    // calc P(P_i = p) and P(P_{~i} = p)
    // time complexity : O(N * I)
    if (calc_fptree_probs(iface_fp_data, tree_bitmasks, in_fptree_prob) < 0)
        return -1;
    // if (calc_fptree_probs(iface_fp_data, tree_bitmask, in_fptree_prob, true) < 0)
    //     return -1;

    // calc P(L_i = l | P_i = p) and P(L_{~i} = l | P_{~i} = p)
    // time complexity : O(N^2 * I)
    if (calc_lm_probs(iface_fp_data, tree_bitmasks) < 0) 
        return -1;

    // calc P(iface = i)
    // time complexity : O(N^2 * I)
    for (uint8_t i = 0; i < this->iface_num; i++) {
        calc_iface_prob(i, iface_fp_data, in_fptree_prob, iface_probs[i]);
        // for joint event '& got into this router'
        iface_probs[i][0] *= in_prob;
        iface_probs[i][1] *= in_prob;
    }

    // calc P(T_{out,i} = t)
    // for joint event '& got into this router
    for (uint8_t i = 0; i < this->iface_num; i++) {

        // std::fill(out_fptree_probs[i].begin(), out_fptree_probs[i].end(), 0.0);
        // out_fptree_probs[i][0] = 1.0;
        calc_out_fptree_prob(i, (*iface_fp_data)[0][i], tree_bitmasks, in_fptree_prob, out_fptree_probs);

        std::cout << "Prob::calc_out_fptree_prob() : P(T_{out," << (int) i << "} = t) :" << std::endl;
        for (uint8_t t = 0; t <= this->n; t++) {
            out_fptree_probs[i][t] *= in_prob;
            std::cout << "\tP(T_{out," << (int) i << "} = " << (int) t << ") = " << out_fptree_probs[i][t] << std::endl;
        }
    }

    return 0;
}

/*
 * \brief   prints the distribution of L_{i,t} random variables, for a single 
 *          iface
 *
 * prints a table in the following format:
 *
 * ---------------------------------------------------------------------------------
 * [|L|]                : | 0              | 1 | ... | |R|_{max} | SUM(P(L_{i,t})) |
 * ---------------------------------------------------------------------------------
 * [P(L_{i,0})]         : | P(L_{i,0} = 0) | . | ... |    ...    |      ...        |
 *     ...              : |      ...       | . | ... |    ...    |      ...        |
 * [P(L_{i,|R|_{max}})] : |      ...       | . | ... |    ...    |      ...        |      
 * ---------------------------------------------------------------------------------
 * 
 * \param   iface   iface for which the L_{i,t} table will be printed
 *
 */
void Prob::print_lm_prob(uint8_t iface, bool iface_complement) {

    // decide if we're calculating the iface or 'iface complement' largest match 
    // probabilities
    std::vector<std::vector<std::vector<__float080> > > * lmp;
    if (!iface_complement) lmp = &(this->lm_cond_prob);
    else lmp = &(this->lm_complement_cond_prob);

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->n + 2); i++)
        printf("-------------");

    printf("\n[|L|]         : |");

    for (uint8_t i = 0; i < (this->n + 1); i++)
        printf(" %-11d|", i);

    if (!iface_complement) printf(" SUM P(L_{i|t}) |");
    else printf(" SUM P(L_{~i|t}) |");

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->n + 2); i++)
        printf("-------------");

    __float080 cumulative_prob = 0.0;
    __float080 prob = 0.0;

    for (uint8_t fptree_size = 0; fptree_size < (this->n + 1); fptree_size++) {

        if (!iface_complement) printf("\n[P(L_{%d|%d}) ] : |", iface, fptree_size);
        else printf("\n[P(L_{~%d|%d})] : |", iface, fptree_size);

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