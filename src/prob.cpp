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
    __float080 lt_ratio) {

    // as a precaution, reset log_prob_fp_neq to 0.0
    std::fill(log_prob_fp_neq.begin(), log_prob_fp_neq.end(), 0.0);

    for (uint8_t f = 0; f < log_prob_fp_neq.size(); f++) {

        // total nr. of entries w/ size f
        __float080 num_entries = iface_fp_data->num_entries * iface_fp_data->f_distr[f] * lt_ratio;
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
    Prob::fp_data * iface_fp_data,
    std::vector<uint8_t> * tree_bitmask) {

    int count = 0;
    for (unsigned int k = 0; k < (*tree_bitmask).size(); k++)
        count += count_set_bits((int) (((*tree_bitmask)[k]) & (iface_fp_data->tree_bitmask[k])));

    return count;
}

__float080 local_tree_ratio(
    Prob::fp_data * iface_fp_data,
    std::vector<uint8_t> * tree_bitmask) {

    // 1) count trees in iface
    int trees = 0, shared_trees = 0;
    for (unsigned int k = 0; k < (*tree_bitmask).size(); k++) {
        trees += count_set_bits((int) iface_fp_data->tree_bitmask[k]);
        shared_trees += count_set_bits((int) (((*tree_bitmask)[k]) & (iface_fp_data->tree_bitmask[k])));
    }

    __float080 lt_ratio = (__float080) (trees - shared_trees) / (__float080) (trees);
    std::cout << "\t [# trees] : " << trees
        << "\n\t [# shared trees] : " << shared_trees
        << "\n\t [lt ratio] : " << lt_ratio << std::endl;

    return lt_ratio;
}

int Prob::calc_lm_prob(
    uint8_t i, 
    Prob::fp_data * iface_fp_data,
    std::vector<uint8_t> * tree_bitmask,
    std::vector<std::vector<__float080> > & out_fptree_probs,
    bool iface_complement) {

    // decide if we're calculating the iface or 'iface complement' largest match 
    // probabilities
    std::vector<std::vector<std::vector<__float080> > > * lmp;
    if (!iface_complement) lmp = &(this->lm_cond_prob);
    else lmp = &(this->lm_complement_cond_prob);

    // aux. arrays used in fp rate calculation:
    //  - P(|FP| =/= f) = (1.0 - P(|FP| = f)) : prob of not having an fp match of size f
    std::vector<__float080> log_prob_fp_neq(n, 0.0);
    //  - P(|FP| <= f) : prob of fp match smaller than or equal to f
    std::vector<__float080> log_prob_fp_smeq(n + 1, 0.0);

    // if there are no entries for this iface, a match is impossible for 
    // this [i][ptree_size][f] combination
    // as such, match size 0 has prob 1.0 of happening
    if (iface_fp_data->num_entries < 1.0) {

        for (int fptree_size = (int) this->n; fptree_size >= 0; fptree_size--) {
            // std::cout << "Prob::calc_lm_prob(iface) : [INFO]: " 
            //     << "\t\n [fptree_size] : " << ((int) fptree_size) << std::endl;
            (*lmp)[i][fptree_size][0] = 1.0;
        }

        return 0;
    }

    // calc the ratio of trees which can origin FP trees locally at the router
    std::cout << "Prob::calc_lm_prob(iface) : [INFO] lt ratio iface  " << (iface_complement ? "~" : "") << (int) i << " :" << std::endl;
    __float080 lt_ratio = local_tree_ratio(iface_fp_data, tree_bitmask);

    std::cout << "Prob::calc_lm_prob(iface) : [INFO] iface " << (iface_complement ? "~" : "") << ((int) i) << " stats :"
        << "\n\t [n, m, k] : " << this->n << ", " << this->m << ", " << this->k
        << "\n\t [# entries] : " << iface_fp_data->num_entries << " (" << (iface_fp_data->num_entries * lt_ratio) << ")";
    //     << "\t\n [lmp size] : " << (*lmp).size()
    //     << "\t\n [f_distr] : ";
    // for(uint8_t f = 0; f < this->n; f++) 
    //     std::cout << "[" << (int) f << "] = " << iface_fp_data->f_distr[f] << ", ";
    std::cout << std::endl;

    // fill P(|FP| =/= f) = (1.0 - P(|FP| = f)), for all f
    calc_log_prob_fp_neq(iface_fp_data, log_prob_fp_neq, lt_ratio);
    // calc P(|FP| <= f), for all f
    for (int f = (int) (this->n - 1); f >= 0; f--)
        log_prob_fp_smeq[f] = log_prob_fp_smeq[f + 1] + log_prob_fp_neq[f];

    if (!iface_complement)
        calc_out_fptree_prob(i, iface_fp_data, log_prob_fp_neq, log_prob_fp_smeq, out_fptree_probs);

    for (int fptree_size = (int) this->n; fptree_size >= 0; fptree_size--) {
        for (int f = (int) this->n; f >= 0; f--) {

            // if f is smaller than the tp size or fp tree size, then f *WON'T 
            // BE* the largest match for sure, so P(L_i = f) = 0.0
            if ((f < iface_fp_data->tp_size) || f < (fptree_size)) {

                (*lmp)[i][fptree_size][f] = 0.0;

            } else {

                // 2 sub-cases (deterministic):
                //  1) f == tp size
                //  2) f > tp size
                if (f == iface_fp_data->tp_size) {

                    // 1) f == tp size
                    // the prob of this event can be calculated as:
                    //  * not having fp matches larger than f : log_prob_fp_smeq[f]
                    (*lmp)[i][fptree_size][f] = exp(log_prob_fp_smeq[f]);

                } else {

                    // 2) f > tp size
                    if (f == fptree_size) {

                        // the prob of this event can be calculated as:
                        //  * not having fp matches larger than f : log_prob_fp_smeq[f]
                        (*lmp)[i][fptree_size][f] = exp(log_prob_fp_smeq[f]);

                    } else {

                        // if the iface is potentially on ptree, then it can't have 
                        // a match larger than the prefix tree. if 
                        // the iface is def. not on ptree, it can be the start of 
                        // a new prefix tree, and as such fp matches 
                        // larger than f are possible.
                        if (iface_fp_data->on_fptree && fptree_size > 0) {

                            (*lmp)[i][fptree_size][f] = 0.0;
                        
                        } else {

                            // the prob of this event can be calculated as:
                            //  * not having fps larger than f : log_prob_fp_smeq[f]
                            //  * AND having a fp of size f : 1.0 - exp(log_prob_fp_neq[f - 1])
                            (*lmp)[i][fptree_size][f] = 
                                exp(log_prob_fp_smeq[f]) 
                                * (1.0 - exp(log_prob_fp_neq[f - 1])) 
                                * ((iface_fp_data->num_entries > 0) ? 1.0 : 0.0);
                        }
                    }
                }
            }
        }
    }

    return 0;
}

int Prob::calc_lm_probs(
    std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
    std::vector<uint8_t> * tree_bitmask,
    std::vector<std::vector<__float080> > & out_fptree_probs) {

    // if (this->has_lm_prob) {
    //     std::cout << "Prob::calc_lm_probs() : [INFO] pre-calculated lm_cond_prob table available. skipping." << std::endl;
    //     return 0;
    // }

    for (uint8_t i = 0; i < this->iface_num; i++) {

        // time complexity : O(N^2 * I)
        calc_lm_prob(i, (*iface_fp_data)[0][i], tree_bitmask, out_fptree_probs);       // iface i
        // time complexity : O(N^2 * I)
        calc_lm_prob(i, (*iface_fp_data)[1][i], tree_bitmask, out_fptree_probs, true); // iface ~i

        print_lm_prob(i);
        print_lm_prob(i, true);
    }

    // mark lm probs as calculated
    this->has_lm_prob = true;

    return 0;
}

// objective : calc P(T_{out,i} = t), for each iface i and all fp tree sizes t
//
// the objective is calculating the probability of the event 'the largest fp tree 
// size a request going out over iface i belongs to equals t', or P(T_{out,i} = t) in 
// short.
// 
// algorithm : 
//  - P(T_{out,i} = t) is influenced by 2 events:
//      1) a fp tree of size t - to which the request is stuck - continues on iface i
//      2) the request falls on a 'fresh' fp tree of size t - due to a local FP 
//         match - which starts on iface i.
//         note: this only happens if the request didn't get stuck to a fp tree on 
//         previous routers. 
//
// time complexity : 
//  - O(N) : 
//      - for each iface i (this function)
//  - O(N * I) : 
//      - for all ifaces
int Prob::calc_out_fptree_prob(
    uint8_t i,
    Prob::fp_data * iface_fp_data,
    std::vector<__float080> log_prob_fp_neq,
    std::vector<__float080> log_prob_fp_smeq,
    std::vector<std::vector<__float080> > & out_fptree_probs) {

    // init all probs to 0.0
    for (uint8_t t = 0; t <= this->n; t++)
        out_fptree_probs[i][t] = 0.0;

    // P(P_{out,i} > 0)
    for (uint8_t t = 1; t <= this->n; t++) {

        // contribution of event 1
        out_fptree_probs[i][t] = this->fptree_prob[0][i][t];
        // contribution of event 2 : local FP match
        // FIXME : isn't '* fptree_prob[0][i][0]' missing here? 
        // i.e the multiplication by P(T_i = 0)...
        out_fptree_probs[i][t] +=
            exp(log_prob_fp_smeq[t])                // not have a FP > t
            * (1.0 - exp(log_prob_fp_neq[t - 1]));  // have a FP of size t
    }

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
    std::vector<uint8_t> * tree_bitmask,
    std::vector<__float080> * in_fptree_prob) {
    // bool iface_complement) {

    // // if iface complement == true, then use the ~i values
    // int ic = 0;
    // if (iface_complement) ic = 1;

    // (1) calculate e[i][t] and tr[i]
    // totals
    __float080 e_total = 0.0;
    __float080 tr_total = 0.0;
    // individual
    std::vector<std::vector<__float080> >  e((this->iface_num), std::vector<__float080>((this->n + 1), 0.0));
    std::vector<__float080> tr(this->iface_num, 0.0);

    for (uint8_t i = 0; i < this->iface_num; i++) {
        // don't consider blocked ifaces
        if ((*iface_fp_data)[0][i]->is_blocked) continue;
        // count shared trees w/ bitmask on iface i
        tr[i] = (__float080) count_shared_trees((*iface_fp_data)[0][i], tree_bitmask);
        // keep track of total # of shared trees
        tr_total += tr[i];

        // % entries in iface i w/ size t
        for (uint8_t t = 1; t < (this->n + 1); t++) {
            e[i][t] = (*iface_fp_data)[0][i]->f_distr[t - 1];
            e_total += e[i][t];
        }
    }

    std::cout << "Prob::calc_fptree_probs() : shared tree ratios :" << std::endl;
    for (uint8_t i = 0; i < this->iface_num; i++)
        std::cout << "\ttr[" << (int) i << "] = " << tr[i] / ((tr_total > 0.0) ? tr_total : 1.0) << std::endl;

    // std::cout << "Prob::calc_fptree_probs() : entry props :" << std::endl;
    // for (uint8_t i = 0; i < this->iface_num; i++)
    //     for (uint8_t t = 0; t < (this->n + 1); t++)
    //         std::cout << "\te[" << (int) i << "][" << (int) t << "] = " << e[i][t] << std::endl;

    // (2) calc P(T_i = t), for all {i,t} pairs, using weight method
    for (uint8_t i = 0; i < this->iface_num; i++) {
        
        // if iface i is blocked, don't bother going on...
        if ((*iface_fp_data)[0][i]->is_blocked) continue;

        // calc P(T_i = t), for all t > 0
        __float080 cumulative_fptree_prob[2] = {0.0, 0.0};
        for (uint8_t t = 1; t < (this->n + 1); t++) {

            this->fptree_prob[0][i][t] = 
                (e[i][t]) 
                * (tr[i] / ((tr_total > 0.0) ? tr_total : 1.0)) 
                * (*in_fptree_prob)[t];

            this->fptree_prob[1][i][t] = ((*in_fptree_prob)[t] - this->fptree_prob[0][i][t]);

            cumulative_fptree_prob[0] += this->fptree_prob[0][i][t];
            cumulative_fptree_prob[1] += this->fptree_prob[1][i][t];
        }

        // set P(T_i = 0) = 1.0 (adjustments are handled in P(L_i = l) calculation))
        this->fptree_prob[0][i][0] = 1.0;
        this->fptree_prob[1][i][0] = 1.0;

        if ((*iface_fp_data)[0][i]->tp_size > 0) {
            this->fptree_prob[1][i][0] = 1.0 - cumulative_fptree_prob[1];
        }

        if ((*iface_fp_data)[1][i]->tp_size > 0) {
            this->fptree_prob[0][i][0] = 1.0 - cumulative_fptree_prob[0];
        }
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

// objective: calc avg nr. of events of type e, for all events
//
// algorithm: 
//  - events are:
//      - EVENT_LLM : choosing iface 0 (local)
//      - EVENT_NLM : prob of not choosing any iface at all, i.e. P(L_i = 0, L_j = 0), for all i, j
//      - EVENT_SLM : sum of P(L_i > L_{~i}), for all i.
//                    P(L_i > L_{~i}) is given by iface_probs[i][0]
//                    time complexity is O(N * I)
//      - EVENT_MLM : 'max. possible nr. events' - P(EVENT_NLM) + P(EVENT_SLM)
//
// time complexity: 
//  - O(I)
//
int Prob::calc_event_num(
    std::vector<std::vector<Prob::fp_data *> > * iface_fp_data,
    std::vector<std::vector<__float080> > iface_probs,
    std::vector<__float080> & event_num) {

    // LLM avg. num_entries: only iface 0 (0 is always the local)
    // FIXME : in this case, we give priority to iface 0, and always consider 
    // P(L_0 >= L_{~0}) probability (instead the 'strict' P(L_0 > L_{~0}))
    event_num[EVENT_LLM] = iface_probs[0][1];
    // NLM & SLM avg. nrs:
    //  - NLM : corresponds to P(L_i = 0 AND L_j = 0, ..., L_n = 0) for all ifaces i, j, ... n
    //          NOTE : this is simply P(L_i = 0 AND L_{~i} = 0), for any i (?)
    event_num[EVENT_NLM] = 0.0;
    //  - SLM : sum P(L_i > L_{~i}), for all i 
    //          these are independent events, i.e. P(L_i > L_{~i}) AND P(L_j > L_{~j}) = 0 (for i != j)
    //          P(L_i > L_{~i}) is given by iface_probs[i][0]
    int blocked_ifaces_num = 0;
    __float080 event_num_total = 0.0;
    for (unsigned int i = 0; i < iface_probs.size(); i++) {
        // count the nr. of block ifaces
        if ((*iface_fp_data)[0][i]->is_blocked) { blocked_ifaces_num++; continue; }
        // (SLM : as above)
        event_num[EVENT_SLM] += iface_probs[i][0];
        // P(L_i >= L_{~i}), will be useful for MLM avg. nr. computation
        event_num_total += iface_probs[i][1];
    }

    // MLM avg. nrs.:
    // FIXME : this still bugs me a lot
    event_num[EVENT_MLM] = event_num_total - (event_num[EVENT_NLM] + event_num[EVENT_SLM]);


    std::cout << "Prob::calc_event_num() : E[event type] = " << std::endl; 
    std::cout << "\tE[EVENT_LLM] = " << event_num[EVENT_LLM] << std::endl;
    std::cout << "\tE[EVENT_NLM] = " << event_num[EVENT_NLM] << std::endl;
    std::cout << "\tE[EVENT_SLM] = " << event_num[EVENT_SLM] << std::endl;
    std::cout << "\tE[EVENT_MLM] = " << event_num[EVENT_MLM] << std::endl;

    return 0;
}

__float080 calc_joint_prob(
    __float080 cond_prob,
    __float080 fptree_prob,
    bool condition) {

    // FIXME: this part is ugly and dubious...
    __float080 _fptree_prob = 0.0;
    if (condition)
        _fptree_prob = 0.0;
    else
        _fptree_prob = fptree_prob;

    return (cond_prob * _fptree_prob);
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
    std::vector<__float080> & iface_prob) {

    //  - P(L_i = l)
    std::vector<std::vector<__float080> > lm_marg_prob((this->iface_num), std::vector<__float080> ((this->n + 1), 0.0));
    //  - P(L_{~i} = l)
    std::vector<std::vector<__float080> > lm_complement_marg_prob((this->iface_num), std::vector<__float080> ((this->n + 1), 0.0));

    // 1) calc P(L_i = l) and P(L_{~i} = l), for all l
    // FIXME : this seems unnecessarily complex...
    for (int f = this->n; f >= 0; f--) {
        for (uint8_t t = 0; t < (this->n + 1); t++) {

            // if iface i has a tp of size s larger or equal to the fp tree 
            // size t (i.e. s >= t), then we don't consider the contribution of 
            // fp matches for P(L_i = s).

            // why? ...

            // calc joint prob P(L_i = l, T_i = t)
            __float080 joint_prob = 
                calc_joint_prob(
                    this->lm_cond_prob[i][t][f],    // P(L_i = l | T_i = t)
                    this->fptree_prob[0][i][t],     // P(T_i = t)
                    (t > 0 && ((*iface_fp_data)[0][i]->tp_size >= t))); // condition
            // add the joint P(L_i = l, T_i = t) to marginal P(L_i = l)
            lm_marg_prob[i][f] += joint_prob;
            // print if P(L_i = l, T_i = t) > 0.0
            if (joint_prob > 0.0) {
                std::cout << "\t[" << (int) t << "] +" 
                    << (joint_prob) << " (" << this->lm_cond_prob[i][t][f] << " x " << this->fptree_prob[0][i][t] << ")"
                    << " to P(L_" << (int) i << " = " << (int) f << ") = " << lm_marg_prob[i][f] << std::endl;
            }

            // calc joint prob P(L_{~i} = l, T_{~i} = t)
            joint_prob = 
                calc_joint_prob(
                    this->lm_complement_cond_prob[i][t][f],     // P(L_{~i} = l | T_{~i} = t)
                    this->fptree_prob[1][i][t],                 // P(T_{~i} = t)
                    (t > 0 && ((*iface_fp_data)[1][i]->tp_size >= t))); // condition
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
    std::vector<uint8_t> * tree_bitmask,    // tree bitmask from prev router
    __float080 in_prob,
    std::vector<__float080> * in_fptree_prob,
    std::vector<std::vector<__float080> > & iface_probs,
    std::vector<__float080> & event_num,
    std::vector<std::vector<__float080> > & out_fptree_probs) {

    // calc P(P_i = p) and P(P_{~i} = p)
    // time complexity : O(N * I)
    if (calc_fptree_probs(iface_fp_data, tree_bitmask, in_fptree_prob) < 0)
        return -1;
    // if (calc_fptree_probs(iface_fp_data, tree_bitmask, in_fptree_prob, true) < 0)
    //     return -1;

    // calc P(L_i = l | P_i = p) and P(L_{~i} = l | P_{~i} = p)
    // time complexity : O(N^2 * I)
    if (calc_lm_probs(iface_fp_data, tree_bitmask, out_fptree_probs) < 0) 
        return -1;

    // calc P(iface = i)
    // time complexity : O(N^2 * I)
    for (uint8_t i = 0; i < this->iface_num; i++) {
        calc_iface_prob(i, iface_fp_data, iface_probs[i]);
        // for joint event '& got into this router'
        iface_probs[i][0] *= in_prob;
        iface_probs[i][1] *= in_prob;
    }

    // calc avg. nr. of events (LLM, NLM, SLM and MLM)
    // time complexity : O(I)
    if (calc_event_num(iface_fp_data, iface_probs, event_num) < 0)
        return -1;

    // for joint event '& got into this router'
    for (uint8_t i = 0; i < this->iface_num; i++) {
        std::cout << "Prob::calc_out_fptree_prob() : P(T_{out," << (int) i << "} = t) :" << std::endl;
        for (uint8_t t = 1; t <= this->n; t++) {
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

    if (!iface_complement) printf(" SUM P(L_{i,t}) |");
    else printf(" SUM P(L_{~i,t}) |");

    printf("\n---------------------");
    for (uint8_t i = 0; i < (this->n + 2); i++)
        printf("-------------");

    __float080 cumulative_prob = 0.0;
    __float080 prob = 0.0;

    for (uint8_t fptree_size = 0; fptree_size < (this->n + 1); fptree_size++) {

        if (!iface_complement) printf("\n[P(L_{%d,%d}) ] : |", iface, fptree_size);
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