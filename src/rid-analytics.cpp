#include <math.h>
#include <cfloat>
#include <algorithm>
#include <string>

#include "tree/tree.hh"
#include "node.h"
#include "graph.h"
#include "rid-analytics.h"

const char * PENALTY_TYPE_STR[] = {"feedback", "fallback"};

double fp_rate(double m, double n, double k, double c) {

    return pow((1.0 - exp(-((n / m) * k))), k * c);
}

int erase_decision_tree(tree<Node *> tree_obj) {

    tree<Node *>::post_order_iterator post_itr = tree_obj.begin_post();
    tree<Node *>::post_order_iterator end_post_itr = tree_obj.end_post();

    while (post_itr != end_post_itr) {

        // delete the Node object associated with the tree node (not the 
        // node itself)
        delete (*post_itr);

        ++post_itr;
    }

    return 0;
}

int is_all_eop(tree<Node *> tree_obj) {

    tree<Node *>::leaf_iterator leaf_itr = tree_obj.begin_leaf();
    tree<Node *>::leaf_iterator end_leaf_itr = tree_obj.end_leaf();

    while (leaf_itr != end_leaf_itr) {

        if ((*leaf_itr)->get_next_tier() != END_OF_PATH) {

            return 0;
        }

        ++leaf_itr;
    }

    return 1;
}

double latency_through_tier(int tier, double * latencies) {

    // you always have to go up to tier 0 and go down (at the end)
    double base_latency = 2.0 * latencies[0];

    // for every additional tier, you do 'valley free routing': 
    // up -> up -> ... -> peer -> down -> down -> ... 
    double penalty_candidate = 0.0;

    // the 'up' & 'down' part
    for (int t = 0; t < tier; t++) {
        penalty_candidate += 2.0 * latencies[t];
    }

    // the 'peer' link latency part
    // FIXME: this isn't always correct but let it stay like that for now
    penalty_candidate += latencies[tier];

    return penalty_candidate + base_latency;
}

double latency_to_content(int * content_sources, double * latencies) {

    double latency_to_content = 2.0 * latencies[0];

    for (int l = 0; l < MAX_TIER_DEPTH; l++) {

         if (content_sources[l] > 0) {
            
            latency_to_content += latencies[l];

            // 'that's it, the rebels are there!'
            break;
         }

         latency_to_content += 2 * latencies[l];
    }

    return latency_to_content;
}

double in_flight_penalty(
    struct RIDAnalytics::rid_analytics_inputs input_params,
    int step,
    int curr_node_tier,
    int curr_node_max_tier,
    double prev_node_latency) {

    double penalty = 0.0;

    if (input_params.fp_resolution_tech == PENALTY_TYPE_FEEDBACK) {

        // with PENALTY_TYPE_FEEDBACK, we first go back to the source of 
        // request...
        penalty = prev_node_latency;

        // ... and then go up to the origin server!
        //penalty += latency_to_content(input_params.content_sources, input_params.latencies);
        penalty += latency_through_tier(input_params.origin_tier, input_params.latencies);

    } else if (input_params.fp_resolution_tech == PENALTY_TYPE_FALLBACK) {

        // with PENALTY_TYPE_FALLBACK, we don't go back to the source of the 
        // request: we immediately forward the request to the closest origin.

        // we use the value calculated in 
        // input_params.penalties[curr_node_max_tier]. however, remember 
        // this is the average distance FROM A WRONG DESTINATION to a correct 
        // source, GIVEN THAT the wrong destination has been reached via 
        // curr_node_max_tier. if we're not at a destination but in the middle 
        // of the way, we need to subtract something...

        // ... for now we determine how far from a network edge we are and 
        // subtract this distance from input_params.penalties[curr_node_max_tier]
        double latency_to_closest_edge = 0.0;

        for (int i = curr_node_tier; i >= 0; i--) {

            latency_to_closest_edge += input_params.latencies[i];
        }

        // HACK: if this is the first time we're at curr_node_tier, and 
        // content becomes visible at curr_node_tier, add 1 unit of 
        // latencies[curr_node_tier] to the result...

        if (input_params.content_sources[curr_node_tier] > 0 && 
            ((step) == curr_node_tier)) {

            // printf("[STEP: %d / @TIER: %d] CS: %d\n", 
            //     step, curr_node_tier, 
            //     input_params.content_sources[curr_node_tier]);
            penalty += input_params.latencies[curr_node_tier];
        }

        penalty += input_params.penalties[curr_node_max_tier] - latency_to_closest_edge;
    }

    return penalty;
}

void calc_penalties(struct RIDAnalytics::rid_analytics_inputs & input_params) {

    input_params.penalties = (double *) calloc(input_params.tier_depth, sizeof(double));

    // furthest tier at which content is available
    for (int i = (input_params.tier_depth - 1); i >= 0; i--) {

        if (input_params.content_sources[i] > 0) {
            
            input_params.origin_tier = i;
            break;
        }
    }

    // get cache tier
    for (int i = 0; i < input_params.tier_depth; i++) {

        if (input_params.content_sources[i] > 0) {

            input_params.cache_tier = i;
            break;
        }
    }

    if (input_params.fp_resolution_tech == PENALTY_TYPE_FEEDBACK) {

        // this case is straightforward: if we're at a wrong destination after 
        // passing through tier m, the penalty will be composed by the 
        // following components:
        //
        //  1) latency of going back to the source of the request, i.e.:
        //      2 * t[0] + 2 * t[1] + ... + 2 * t[m - 1] + t[m] 
        //      
        //      note we go through each tier i < m twice (path from wrong 
        //      destination @layer 0 to layer m + path from layer m to 
        //      destination), hence the 2 * t[i < m] 
        //
        //  2) latency of going up to the first [v]isible origin server which 
        //      is able to serve the content:
        //      2 * t[0] + 2 * t[1] + ... + 2 * t[v - 1] + t[v]
        //

        // calculate a penalty for each tier max. m
        for (int m = 0; m < input_params.tier_depth; m++) {

            // 1) going back through tier m to content source
            input_params.penalties[m] = latency_through_tier(m, input_params.latencies);

            // 2) we got to the source of the request. now we move up, till 
            // a content source is visible
            input_params.penalties[m] += latency_through_tier(input_params.origin_tier, input_params.latencies);
        }

    } else if (input_params.fp_resolution_tech == PENALTY_TYPE_FALLBACK) {

        double amortized_penalty = 0.0;
        double prob_prev = 1.0, prob_next = 1.0;

        for (int m = 0; m < input_params.tier_depth; m++) {

            // here's the logic: if went up as high as tier m to get 
            // to a wrong destination AND there are visible sources at tier m, 
            // then one of the following might be true:
            //  a) the correct server might be in the same domain in 
            //      tier i < m domain as the wrong server
            //  b) the correct server might be in a different domain at tier 
            //      m
            //
            // the penalty of case a) is smaller than b)
            if (input_params.content_sources[m] > 0) {

                // the correct server is accessible via tier i, and so 
                // penalty should be latency_through_tier(i, latencies)

                // there is a probability associated with each penalty. how do 
                // we determine it?
                //  -# each tier i has domains[i] domains (on avg. provided as input)
                //  -# there are content_sources[m], which can be thought as 
                //      uniformly distributed (ASSUMPTION) over all domains 

                prob_prev = 1.0;
                prob_next = 1.0;

                amortized_penalty = 0.0;

                for (int i = m; i >= 0; i--) {
                    
                    // probabilities always < 1.0 please: if there are more 
                    // content sources than domains (and assuming a uniform 
                    // distribution of sources per domains), we truncate 
                    // probabilities here...
                    prob_next = min(((double) input_params.content_sources[m]) / ((double) (input_params.domains[i])), 1.0);
                    amortized_penalty += prob_prev * (1.0 - prob_next) * latency_through_tier(i, input_params.latencies);

                    // printf("penalties[%d.%d] = %f * %f = %f\n", 
                    //     m, i,
                    //     (1.0 - prob_next), 
                    //     latency_through_tier(i, input_params.latencies), 
                    //     amortized_penalty);

                    prob_prev = prob_prev * prob_next;
                }

                // special operation for the minimum latency case...
                // FIXME: not sure if totally correct
                amortized_penalty += prob_prev * 3.0 * input_params.latencies[0];

                // printf("penalties[%d] += %f * 3.0 * %f = %f\n", 
                //     m,
                //     prob_prev, 
                //     input_params.latencies[0], 
                //     amortized_penalty);

            } else {

                amortized_penalty = latency_to_content(
                        input_params.content_sources,
                        input_params.latencies);
            }

            input_params.penalties[m] = amortized_penalty;
        }
    }
}

void run_checksum(
    tree<Node *> decision_tree, 
    struct RIDAnalytics::rid_analytics_inputs input_params) {

    tree<Node *>::leaf_iterator leaf_itr = decision_tree.begin_leaf();
    tree<Node *>::leaf_iterator end_leaf_itr = decision_tree.end_leaf();

    printf("\n*** INTERMEDIATE CHECKSUM (%s | %s | %d.%d) ***\n\n", 
        input_params.title, 
        PENALTY_TYPE_STR[input_params.fp_resolution_tech],
        input_params.cache_tier, input_params.origin_tier);

    double avg_latency = 0.0;
    double checksum = 0.0;
    char * node_str = NULL;

    while (leaf_itr != end_leaf_itr) {

        if ((*leaf_itr)->get_prob_val() > 0.0) {

            node_str = (*leaf_itr)->to_string();

            printf("[DEPTH %d | @TIER %d] : %s\n", 
                decision_tree.depth(leaf_itr), 
                (*leaf_itr)->get_curr_tier(),
                node_str);

            free(node_str);

            checksum += (*leaf_itr)->get_prob_val();
            avg_latency += ((*leaf_itr)->get_prob_val() * (*leaf_itr)->get_latency_val());
        }

        ++leaf_itr;
    }

    printf("\n[CHECKSUM : %-.8E]\n", checksum);
    printf("[AVG_LATENCY : %-.8E]\n", avg_latency);

    // wait for user input to continue...
    printf("\npress 'c' to continue\n");

    char entered = 0;
    while (entered != 'c') { 
        entered = getchar(); 
    }
}

void print_params(struct RIDAnalytics::rid_analytics_inputs input_params) {

    printf("\n*** PARAMETERS ***\n");

    // initial row w/ tier numbers
    printf("\n\n[TIER]    : |");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf(" %-8d|", i);

    // line separating header row and table body
    printf("\n-------------");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf("----------");

    // fp probability per tier
    printf("\n[FP-PROB] : |");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf(" %-.2E|", input_params.fp_prob[i]);

    // latencies to each tier
    printf("\n[LATENCY] : |");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf(" %-.6f|", input_params.latencies[i]);

    // penalties per tier
    printf("\n[PENALTY] : |");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf(" %-.6f|", input_params.penalties[i]);

    // alpha values at each tier
    printf("\n[ALPHA]   : |");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf(" %-.2E|", input_params.alpha[i]);

    printf("\n\n");
}

int RIDAnalytics::run_model(
    double & avg_latency,
    struct RIDAnalytics::rid_analytics_inputs input_params,
    unsigned int modes, 
    FILE ** outcomes_file,
    int outcomes_files_size,
    char * title,
    std::string data_dir) {

    // // latency values per tier, in multiples of some time unit.
    // double latencies[MAX_TIER_DEPTH];
    // // FIXME: just initialize everything to 1.0
    // for (int tier = 0; tier < input_params.tier_depth; tier++)
    //     latencies[tier] = 1.0;

    // // point the latencies array in the input params struct to the array 
    // // just created
    // input_params.latencies = latencies;

    // the penalty array. our goal is to calculate the relaying 
    // penalty - PENALTY(m) - for a request, given the following conditions: 
    //
    //  -# the request has crossed a tier as high as tier m
    //  -# we know the distribution of sources per tier for the scenario 
    //     being evaluated
    //  -# we know the avg. latencies of forwarding between routers at tier l 
    //
    // this can be pre-calculated 'a priori'. how? see below...
    //
    calc_penalties(input_params);

    // print all input parameters (if verbose mode is set)
    if ((modes & MODE_VERBOSE)) {
        print_params(input_params);
    }

    // ************************************************************************
    // build the n-ary tree representing all possible forwarding decisions 
    // and associated outcomes in the path of a request towards a content 
    // source, for the scenario(s) described in the .scn file(s)
    // ************************************************************************
    tree<Node *> decision_tree;
    tree<Node *>::iterator root;

    // initialize the probability tree using begin() and add a root 
    // (dummy) node with insert(). after this, we use append_child() all the 
    // way.
    Node * root_node = new Node(0, 0, 0, 1.0, input_params.latencies[0], Node::ORI_NODE);

    root = decision_tree.begin();
    root = decision_tree.insert(root, root_node);

    // we fill the tree tier by tier. we add children to each non-leaf node 
    // of the previous tier. to iterate over nodes of a given tier, 
    // we use fixed_depth_iterators of tree.hh.
    tree<Node *>::fixed_depth_iterator depth_itr;
    tree<Node *>::fixed_depth_iterator end_depth_itr;

    Node * mhs_node = NULL;
    Node * mhd_node = NULL;
    Node * fpo_node = NULL;
    Node * tpo_node = NULL;
    Node * def_node = NULL;

    int curr_node_tier = 0;
    int curr_node_max_tier = 0;
    double prev_node_prob = 1.0;
    double prev_node_latency = 0.0;
    Node::Type prev_node_type = Node::UNKNOWN;

    double path_latency = 0.0;
    double penalty_latency = 0.0;

    double pos_prob_local = 1.0;
    double alpha_local = 1.0;
    double single_fp_prob = 1.0;

    // start the .dot file for rendering the probability tree
    Graph * graphviz_graph = NULL;

    if ((modes & MODE_SAVEGRAPH)) {

        graphviz_graph = new Graph(std::string(data_dir + "/" + title).c_str());        
    }

    // depth and breadth are sort of (x, y) coordinates for nodes in a tree. 
    // these will be used to generate node names for the .dot file. the origin 
    // of the tree has coordinates (0, 0).
    int depth = 0;
    int breadth = 0;

    // note the limit condition : on the longest possible path (i.e. going 
    // through the max. tier) we make 2 * tier_depth decisions
    for (depth = 0; depth < ((2 * input_params.tier_depth)); depth++) {

        // update the iterators
        depth_itr = decision_tree.begin_fixed(root, depth);

        breadth = 0;

        // cycle through the previous decisions, append children
        while (decision_tree.is_valid(depth_itr)) {

            curr_node_tier = (*depth_itr)->get_next_tier();
            curr_node_max_tier = (*depth_itr)->get_max_tier();

            // if this is a 'dead end', don't append children to it
            if (curr_node_tier == END_OF_PATH) {
                
                // obviously, increment the fixed_depth_iterator so that we 
                // check the next node at the same depth
                ++depth_itr;

                // increment the breadth coordinate
                breadth++;

                continue;
            }

            // save the information about the node to which we will append
            // the children
            prev_node_prob = (*depth_itr)->get_prob_val();
            prev_node_latency = (*depth_itr)->get_latency_val();
            prev_node_type = (*depth_itr)->get_type();

            // create the children Node objects, specifying the diff types
            mhs_node = new Node(curr_node_max_tier, curr_node_tier, Node::MHS_NODE);
            mhd_node = new Node(curr_node_max_tier, curr_node_tier, Node::MHD_NODE);
            tpo_node = new Node(curr_node_max_tier, curr_node_tier, Node::TPO_NODE);
            fpo_node = new Node(curr_node_max_tier, curr_node_tier, Node::FPO_NODE);
            def_node = new Node(curr_node_max_tier, curr_node_tier, Node::DEF_NODE);

            // for each node type, we answer 4 questions, which determine the 
            // values associated with nodes and transitions:
            //
            //  Q1) what's the 'TP reality'? no TPs in the table? 1 single TP? 
            //      n TPs? this will determine the way we handle FP 
            //      probabilities and calculate transition probabilities between 
            //      states
            //  Q2) in what tier will the next node be? going up a tier? down? 
            //      peering? remember we use 'valley free routing' in a 
            //      hierarchy of tiers.
            //  Q3) what's the latency a certain lookup outcome implies? this 
            //      includes the correct calculation of penalties due to 
            //      relaying.
            //  Q4) what's the outcome of the sequence of decisions at the 
            //      end of the path? e.g. is the request packet delivered 
            //      correctly?

            // the MHD node is a special case in our model. in terms of 
            // latency, next tier and outcomes, it is always treated in the 
            // same way, independently of the source state (the probabilities 
            // are dependent on the source though...)
            
            // Q2) for MHD decisions, penalty latencies are always 
            // applied. however, we can't use the penalty array directly, since 
            // the relaying here MAY be made when the packet is in-flight (and 
            // thus at some tier i) and not from some wrong origin server (at 
            // tier 0).
            double mhd_penalty = in_flight_penalty(
                    input_params,
                    depth,
                    curr_node_tier,
                    curr_node_max_tier,
                    prev_node_latency);

            mhd_node->set_latency_val(prev_node_latency + mhd_penalty);

            // printf("[@tier %d / %d] mhd_node->set_latency_val(%f + %f = %f)\n",
            //         curr_node_tier,
            //         curr_node_max_tier, 
            //         prev_node_latency,
            //         mhd_penalty, 
            //         prev_node_latency + mhd_penalty);

            // Q3) MHD decisions always end up in relaying the packet
            // towards an origin server
            mhd_node->set_next_tier(END_OF_PATH);
            mhd_node->set_outcome(OUTCOME_RELAYED);

            if (prev_node_type == Node::ORI_NODE || prev_node_type == Node::DEF_NODE) {     

                // Q1) TP reality

                // NOTE: we use 'curr_node_max_tier' as the index to determine 
                // the 'TP reality'. why?
                // 
                // we cannot use 'curr_node_tier' as the index: the 
                // content_sources array only tells us about the content 
                // sources than become 'visible' upon the first time a packet
                // enters into tier i, i.e. within the domains 
                // at tier i which are closer to the source of a request. when 
                // we're following a 'down' link, the domains at tier i we 
                // encounter are not the same as those we encounter when going 
                // up. nevertheless, these will inherit the presence 
                // of TPs or FPs from upper levels.
                if (input_params.content_sources[curr_node_max_tier] == 0) {

                    // TP reality = no TPs in forwarding table
                    // 
                    //  NEXT STATE | HAPPENS WHEN?                      
                    // ------------|------------------------------------ 
                    //  MHD        | (FPs light up) AND (iface is diff) 
                    //  FPO        | (FPs light up) AND (iface is same)
                    //  MHS        | Never
                    //  DEF        | NOT (FPs light up)

                    // *** probabilities ***

                    // FIXME #1: single_fp_prob is incorrect. 
                    // FIXME #1: besides, taking it into account makes me think 
                    // this particular case doesn't make any sense.
                    
                    // FIXME #2: ok, i've decided to only consider the case of 
                    // only having FPs and not necessarily a single one. this 
                    // makes sense: what we want to check is the fact that 
                    // we only have FPs to choose from which will forcefully 
                    // lead us in the wrong direction.

                    // why is this good? well, it frees me from having to calc
                    // the SFP probability which is a pain to do without 
                    // access to |F\R| and table sizes.
                    //single_fp_prob = 0.001 * input_params.fp_prob[curr_node_tier];
                    single_fp_prob = 0.0;
                    pos_prob_local = (input_params.fp_prob[curr_node_tier] - single_fp_prob);
                    //mhs_node->set_prob_val(prev_node_prob * (pos_prob_local * (input_params.alpha[curr_node_tier])));
                    mhd_node->set_prob_val(prev_node_prob * (pos_prob_local * (1.0 - input_params.alpha[curr_node_tier])));
                    
                    // FIXME #2: i basically replace the probability of an FPO 
                    // outcome with that of a MHS: again, if we have different 
                    // interfaces, we detect (at least) one FP, and thus relay 
                    // the packet
                    //fpo_node->set_prob_val(prev_node_prob * single_fp_prob);
                    fpo_node->set_prob_val(prev_node_prob * (pos_prob_local * (input_params.alpha[curr_node_tier])));
                    def_node->set_prob_val(prev_node_prob * (1.0 - input_params.fp_prob[curr_node_tier]));

                    // IMPOSSIBLE CASES
                    tpo_node->set_prob_val(0.0);
                    mhs_node->set_prob_val(0.0);

                    // *** latencies *** 

                    //mhs_node->set_latency_val(prev_node_latency + input_params.latencies[curr_node_tier]);
                    fpo_node->set_latency_val(prev_node_latency + input_params.latencies[curr_node_tier]);
                    def_node->set_latency_val(prev_node_latency + input_params.latencies[curr_node_tier]);                    

                    // Q3) when coming from a DEF lookup state, packets will 
                    // always be forwarded 'up' or to a 'peer' domain. packets 
                    // can also be relayed to an origin server if the MHD state 
                    // is reached. the rules are the following:

                    //  NEXT STATE      | NEXT TIER                      
                    // -----------------|------------------------------------
                    //  MHD             | relayed to orig. server
                    //  MHS, TPO, FPO   | current tier
                    //  DEF             | current tier + 1 (or dropped)
                    //
                    // the above applies to all TP realities

                    // FIXME #2: MHS changes to END_OF_PATH
                    //mhs_node->set_next_tier(curr_node_tier);
                    mhs_node->set_next_tier(END_OF_PATH);
                    fpo_node->set_next_tier(curr_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    
                    // Q3) default route : go up a tier or (drop the packet if 
                    // we're already at the top and no TPs exist)
                    if (curr_node_max_tier == (input_params.tier_depth - 1)) {

                        def_node->set_next_tier(END_OF_PATH);
                        def_node->set_outcome(OUTCOME_DROPPED);

                    } else {

                        def_node->set_next_tier(curr_node_tier + 1);
                        def_node->set_max_tier(curr_node_tier + 1);
                    }

                    // add nodes to the .dot file for graph rendering
                    if ((modes & MODE_SAVEGRAPH)) {

                        //graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (pos_prob_local * (input_params.alpha[curr_node_tier])));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (pos_prob_local * (1.0 - input_params.alpha[curr_node_tier])));
                        graphviz_graph->add_node((*depth_itr), depth, fpo_node, depth + 1, breadth, (pos_prob_local * (input_params.alpha[curr_node_tier])));
                        graphviz_graph->add_node((*depth_itr), depth, def_node, depth + 1, breadth, (1.0 - input_params.fp_prob[curr_node_tier]));

                        //Node * nodes[] = {mhs_node, mhd_node, fpo_node, def_node};
                        Node * nodes[] = {mhd_node, fpo_node, def_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_max_tier] > 0) {

                    // TP reality = only 1 TP in forwarding table
                    // 
                    //  NEXT STATE | HAPPENS WHEN?                      
                    // ------------|------------------------------------
                    //  MHS        | (FPs light up) AND (iface is same) 
                    //  MHD        | (FPs light up) AND (iface is diff) 
                    //  TPO        | NOT (FPs light up)                 

                    // *** probabilities ***
                    pos_prob_local = input_params.fp_prob[curr_node_tier];
                    mhs_node->set_prob_val(prev_node_prob * (pos_prob_local * (input_params.alpha[curr_node_tier])));
                    mhd_node->set_prob_val(prev_node_prob * (pos_prob_local * (1.0 - input_params.alpha[curr_node_tier])));
                    tpo_node->set_prob_val(prev_node_prob * (1.0 - pos_prob_local));
                    // IMPOSSIBLE CASES
                    fpo_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + input_params.latencies[curr_node_tier]);
                    tpo_node->set_latency_val(prev_node_latency + input_params.latencies[curr_node_tier]);

                    // *** next tier *** 
                    mhs_node->set_next_tier(curr_node_tier);
                    tpo_node->set_next_tier(curr_node_tier);
                    fpo_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (pos_prob_local * (input_params.alpha[curr_node_tier])));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (pos_prob_local * (1.0 - input_params.alpha[curr_node_tier])));
                        graphviz_graph->add_node((*depth_itr), depth, tpo_node, depth + 1, breadth, (1.0 - pos_prob_local));

                        Node * nodes[] = {mhs_node, mhd_node, tpo_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::FPO_NODE) {

                // **********************
                // *** GENERAL VALUES ***
                // **********************

                // Q2) in general, if the previous state isn't DEF, we will now 
                // go down a tier
                int next_node_tier = curr_node_tier - 1;

                // Q2) check if the request is about to be delivered to a 
                // content source
                double penalty_latency_fpo = 0.0;
                double penalty_latency_mhs = 0.0;

                if (next_node_tier == END_OF_PATH) {

                    path_latency = input_params.latencies[curr_node_tier];

                    // Q1) quick-fix for the fp probability value to use in 
                    // the probability calculations (notice that if 
                    // next_node_tier == -1 we would have a segfault)
                    pos_prob_local = input_params.fp_prob[curr_node_tier];
                    alpha_local = input_params.alpha[curr_node_tier];
                    
                    // Q3) the penalty latency of an END_OF_PATH fpo node shall 
                    // be set according to the pre-calculated penalty
                    penalty_latency_fpo = input_params.penalties[curr_node_max_tier];

                    // Q4) the outcome of END_OF_PATH fpo nodes shall ALWAYS be 
                    // OUTCOME_ICACHE (incorrect delivery)
                    fpo_node->set_outcome(std::string(OUTCOME_ICACHE));

                    // Q3 & Q4) if TPs are involved, the latency and outcome of 
                    // a MHS node shall ALWAYS be 0.0 and OUTCOME_CCACHE
                    if (input_params.content_sources[curr_node_max_tier] > 0) {

                        mhs_node->set_outcome(std::string(OUTCOME_CCACHE));
                    
                    } else {

                        penalty_latency_mhs = input_params.penalties[curr_node_max_tier];
                        mhs_node->set_outcome(std::string(OUTCOME_ICACHE));
                    }

                } else {

                    pos_prob_local = input_params.fp_prob[next_node_tier];
                    alpha_local = input_params.alpha[next_node_tier];

                    // the next component of path latency will be that of 
                    // the tier below: notice that we're essentially forwarding 
                    // between hops within the same tier (in this case, 
                    // next_node_tier)
                    path_latency = input_params.latencies[next_node_tier];
                } 

                // Q1) which TP reality are we in ? no TPs? single TP? 
                // multiple TPs? 
                if (input_params.content_sources[curr_node_max_tier] == 0) {

                    // TP reality = no TPs in forwarding table
                    // 
                    //  NEXT STATE | HAPPENS WHEN?                      
                    // ------------|------------------------------------
                    //  MHS        | ((FPs light up) AND (NOT (only 1 FP))) AND (iface is same) 
                    //  MHD        | ((FPs light up) AND (NOT (only 1 FP))) AND (iface is diff) 
                    //  FPO        | (only 1 FP)
                    //  DEF        | NOT (FPs light up)

                    // *** probabilities ***

                    // FIXME #2 (see above)
                    //single_fp_prob = 0.001 * input_params.fp_prob[curr_node_tier];
                    single_fp_prob = 0.0;
                    //mhs_node->set_prob_val(prev_node_prob * ((pos_prob_local - single_fp_prob) * alpha_local));
                    mhd_node->set_prob_val(prev_node_prob * ((pos_prob_local - single_fp_prob) * (1.0 - alpha_local)));
                    fpo_node->set_prob_val(prev_node_prob * ((pos_prob_local - single_fp_prob) * alpha_local));
                    def_node->set_prob_val(prev_node_prob * (1.0 - pos_prob_local));

                    // IMPOSSIBLE CASES
                    tpo_node->set_prob_val(0.0);
                    mhs_node->set_prob_val(0.0);

                    // *** latencies *** 


                    //mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency_mhs);
                    fpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency_fpo);
                    def_node->set_latency_val(prev_node_latency  
                                                + in_flight_penalty(
                                                    input_params,
                                                    depth,
                                                    curr_node_tier,
                                                    curr_node_max_tier,
                                                    prev_node_latency));

                    // *** next tier *** 

                    mhs_node->set_next_tier(END_OF_PATH);
                    fpo_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);

                    // Q4) a packet which is going down and gets into a 
                    // DEF state is dropped
                    def_node->set_next_tier(END_OF_PATH);
                    def_node->set_outcome(OUTCOME_DROPPED);

                    if ((modes & MODE_SAVEGRAPH)) {

                        //graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, ((pos_prob_local - single_fp_prob) * alpha_local));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, ((pos_prob_local - single_fp_prob) * (1.0 - alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, fpo_node, depth + 1, breadth, (pos_prob_local * (alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, def_node, depth + 1, breadth, (1.0 - pos_prob_local));

                        //Node * nodes[] = {mhs_node, mhd_node, fpo_node, def_node};
                        Node * nodes[] = {mhd_node, fpo_node, def_node};                        
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_max_tier] > 0) {

                    // TP reality = 1 or more TPs in RID table
                    // 
                    // this is impossible by definition. if there are TPs we 
                    // cannot only have FPs.
                    //
                    // FIXME: this needs a bit more of thinking... in what 
                    // case could a transition between FPO -> TPO happen? 
                    // caching? 

                    // *** probabilities ***
                    mhs_node->set_prob_val(0.0);
                    mhd_node->set_prob_val(0.0);
                    tpo_node->set_prob_val(0.0);
                    fpo_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** next tier *** 
                    mhs_node->set_next_tier(END_OF_PATH);
                    tpo_node->set_next_tier(END_OF_PATH);
                    fpo_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);
                }

            } else if (prev_node_type == Node::TPO_NODE) {

                // **********************
                // *** GENERAL VALUES ***
                // **********************

                // Q2) in general, if the previous state isn't DEF, we will now 
                // go down a tier
                int next_node_tier = curr_node_tier - 1;

                // Q3) for transitions from TPO to {MHS, TPO}, penalty latencies
                // are 0.0
                penalty_latency = 0.0;

                // Q4) likewise, deliveries are always correct for {MHS, TPO} 
                // transitions for END_OF_PATH
                if (next_node_tier == END_OF_PATH) {

                    path_latency = input_params.latencies[curr_node_tier];

                    mhs_node->set_outcome(std::string(OUTCOME_CCACHE));
                    tpo_node->set_outcome(std::string(OUTCOME_CCACHE));

                    pos_prob_local = input_params.fp_prob[curr_node_tier];
                    alpha_local = input_params.alpha[curr_node_tier];

                } else {

                    // the next component of path latency will be that of 
                    // the tier below: notice that we're essentially forwarding 
                    // between hops within the same tier (in this case, 
                    // next_node_tier)
                    path_latency = input_params.latencies[next_node_tier];

                    pos_prob_local = input_params.fp_prob[next_node_tier];
                    alpha_local = input_params.alpha[next_node_tier];
                } 

                // Q1) TP reality
                if (input_params.content_sources[curr_node_max_tier] == 0) {         

                    // TP reality = no TPs
                    // if we come from TPO state, there must be at least 1 
                    // TP in the table. therefore, transitions within the 
                    // 'no TP' reality are impossibru

                    // *** probabilities ***
                    mhs_node->set_prob_val(0.0);
                    mhd_node->set_prob_val(0.0);
                    tpo_node->set_prob_val(0.0);
                    fpo_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** next tier *** 
                    mhs_node->set_next_tier(END_OF_PATH);
                    mhd_node->set_next_tier(END_OF_PATH);
                    tpo_node->set_next_tier(END_OF_PATH);
                    fpo_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                } else if (input_params.content_sources[curr_node_max_tier] > 0) {   

                    // TP reality = 1 or more TPs in forwarding table
                    // 
                    //  NEXT STATE | HAPPENS WHEN?                      
                    // ------------|------------------------------------
                    //  MHS        | (FPs light up) AND (iface is same) 
                    //  MHD        | (FPs light up) AND (iface is diff) 
                    //  TPO        | NOT (FPs light up)                 

                    // *** probabilities ***
                    // non-zero probabilities
                    mhs_node->set_prob_val(prev_node_prob * (pos_prob_local * alpha_local));
                    mhd_node->set_prob_val(prev_node_prob * (pos_prob_local * (1.0 - alpha_local)));
                    tpo_node->set_prob_val(prev_node_prob * (1.0 - pos_prob_local));
                    // everything else is impossible...
                    fpo_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(next_node_tier);
                    fpo_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (pos_prob_local * alpha_local));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (pos_prob_local * (1.0 - alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, tpo_node, depth + 1, breadth, (1.0 - pos_prob_local));

                        Node * nodes[] = {mhs_node, mhd_node, tpo_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                }

            } else if (prev_node_type == Node::MHS_NODE) {

                // Q2) in general, if the previous state isn't DEF, we will now 
                // go down a tier
                int next_node_tier = curr_node_tier - 1;

                // Q3 & Q4) if the request is to be delivered and TPs are involved, 
                // the delivery is correct and penalty is 0.0
                if (next_node_tier == END_OF_PATH) {

                    // penalty due to relaying: if TPs are involved, no 
                    // latency is added and the packet is delivered correctly
                    path_latency = input_params.latencies[curr_node_tier];

                    if (input_params.content_sources[curr_node_max_tier] > 0) {

                        penalty_latency = 0.0;
                        mhs_node->set_outcome(std::string(OUTCOME_CCACHE));
                        tpo_node->set_outcome(std::string(OUTCOME_CCACHE));

                    } else {
                        
                        // ... otherwise
                        penalty_latency = input_params.penalties[curr_node_max_tier];
                        mhs_node->set_outcome(std::string(OUTCOME_ICACHE));
                        fpo_node->set_outcome(std::string(OUTCOME_ICACHE));
                    }

                    pos_prob_local = input_params.fp_prob[curr_node_tier];
                    alpha_local = input_params.alpha[curr_node_tier];

                } else {

                    // the next component of path latency will be that of 
                    // the tier below: notice that we're essentially forwarding 
                    // between hops within the same tier (in this case, 
                    // next_node_tier)
                    path_latency = input_params.latencies[next_node_tier];
                    penalty_latency = 0.0;

                    pos_prob_local = input_params.fp_prob[next_node_tier];
                    alpha_local = input_params.alpha[next_node_tier];
                } 

                // Q3) which TP reality are we in ? no TPs? single TP? 
                // multiple TPs? 
                if (input_params.content_sources[curr_node_max_tier] == 0) {

                    // FIXME #3: this should be impossible. if 
                    // we now consider MHS as being composed by both TPs and 
                    // FPs (non-exclusively), it is somewhat impossible to 
                    // get to this state. this corresponds to the case where 
                    // we previously had TPs in the table and then suddenly 
                    // they disappear... 

                    // TP reality = no TPs in forwarding table
                    // 
                    //  NEXT STATE | HAPPENS WHEN?                      
                    // ------------|------------------------------------
                    //  FPO        | (FPs light up) AND (iface is same) 
                    //  MHD        | (FPs light up) AND (iface is diff) 
                    //  MHS        | Never
                    //  DEF        | NOT (FPs light up)

                    // *** probabilities ***

                    // FIXME #2: (see above)
                    //single_fp_prob = 0.001 * pos_prob_local;
                    single_fp_prob = 0.0;
                    //mhs_node->set_prob_val(prev_node_prob * ((pos_prob_local - single_fp_prob) * alpha_local));
                    mhd_node->set_prob_val(prev_node_prob * ((pos_prob_local - single_fp_prob) * (1.0 - alpha_local)));
                    fpo_node->set_prob_val(prev_node_prob * ((pos_prob_local - single_fp_prob) * alpha_local));
                    def_node->set_prob_val(prev_node_prob * (1.0 - pos_prob_local));
                    // tpos are impossible (by definition)
                    tpo_node->set_prob_val(0.0);
                    mhs_node->set_prob_val(0.0);

                    // *** latencies *** 

                    //mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    fpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency  
                                                + in_flight_penalty(
                                                    input_params,
                                                    depth,
                                                    curr_node_tier,
                                                    curr_node_max_tier,
                                                    prev_node_latency));                    

                    // *** next tier *** 

                    //mhs_node->set_next_tier(next_node_tier);
                    mhs_node->set_next_tier(END_OF_PATH);
                    fpo_node->set_next_tier(next_node_tier);
                    // Q4) a packet which is going down and gets into a DEF 
                    // state is dropped
                    def_node->set_next_tier(END_OF_PATH);
                    def_node->set_outcome(OUTCOME_DROPPED);
                    tpo_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {

                        //graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, ((pos_prob_local - single_fp_prob) * alpha_local));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, ((pos_prob_local - single_fp_prob) * (1.0 - alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, fpo_node, depth + 1, breadth, (pos_prob_local * (alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, def_node, depth + 1, breadth, (1.0 - pos_prob_local));

                        //Node * nodes[] = {mhs_node, mhd_node, fpo_node, def_node};
                        Node * nodes[] = {mhd_node, fpo_node, def_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_max_tier] > 0) {

                    // TP reality = 1 single TP in forwarding table
                    // 
                    //  NEXT STATE | HAPPENS WHEN?                      
                    // ------------|------------------------------------
                    //  MHS        | (FPs light up) AND (iface is same) 
                    //  MHD        | (FPs light up) AND (iface is diff) 
                    //  TPO        | NOT (FPs light up)

                    // *** probabilities ***
                    mhs_node->set_prob_val(prev_node_prob * (pos_prob_local * (alpha_local)));
                    mhd_node->set_prob_val(prev_node_prob * (pos_prob_local * (1.0 - alpha_local)));
                    tpo_node->set_prob_val(prev_node_prob * (1.0 - pos_prob_local));
                    // impossible events...
                    fpo_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(next_node_tier);
                    fpo_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {

                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (pos_prob_local * (alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (pos_prob_local * (1.0 - alpha_local)));
                        graphviz_graph->add_node((*depth_itr), depth, tpo_node, depth + 1, breadth, (1.0 - pos_prob_local));

                        Node * nodes[] = {mhs_node, mhd_node, tpo_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::MHD_NODE) {

                // EDIT: we should never get in here...
            }

            // finally append the 3 children
            decision_tree.append_child(depth_itr, mhs_node);
            decision_tree.append_child(depth_itr, mhd_node);
            decision_tree.append_child(depth_itr, fpo_node);
            decision_tree.append_child(depth_itr, tpo_node);
            decision_tree.append_child(depth_itr, def_node);

            breadth++;

            // go to the next node in the previous tier and repeat the process
            ++depth_itr;

            // run a sanity check on the nodes just added...
            if ((modes & MODE_VERBOSE))
                run_checksum(decision_tree, input_params);
        }

        if (is_all_eop(decision_tree) > 0)
            break;
    }

    // wrap the .dot file up
    if ((modes & MODE_SAVEGRAPH))
        graphviz_graph->terminate();

    // the decision tree is finished. let's iterate through the leafs to learn 
    // a bit about latencies for this scenario.
    tree<Node *>::leaf_iterator leaf_itr = decision_tree.begin_leaf();
    tree<Node *>::leaf_iterator end_leaf_itr = decision_tree.end_leaf();

    // print some results for quick checking
    if ((modes & MODE_VERBOSE))
        printf("\n*** FINAL RESULTS ***\n\n");

    double cache_latency = DBL_MAX;
    double checksum = 0.0;
    
    char * node_str = NULL;

    while (leaf_itr != end_leaf_itr) {

        if ((*leaf_itr)->get_prob_val() > 0.0) {

            if ((modes & MODE_VERBOSE)) {

                node_str = (*leaf_itr)->to_string();

                printf("[DEPTH %d] : %s\n", 
                    decision_tree.depth(leaf_itr), 
                    node_str);

                free(node_str);
            }

            checksum += (*leaf_itr)->get_prob_val();
            avg_latency += ((*leaf_itr)->get_prob_val() * (*leaf_itr)->get_latency_val());

            if (cache_latency > (*leaf_itr)->get_latency_val())
                cache_latency = (*leaf_itr)->get_latency_val();

            if ((modes & MODE_SAVEOUTCOMES)) {

                // if requested, save outcome information in files pointed 
                // by the outcomes_file array. this will basically save lines 
                // in the format: 
                // [fp@<_fp_tier>],[_fp_prob[<_fp_tier>]],[a@<_alpha_tier>],[_alpha[_alpha_tier]],[penalty_type],[FP RESOLUTION TECH],[LOOKUP OUTCOME TYPE],[OUTCOME TYPE PROBABILITY OF OCCURRENCE]
                for (int i = 0; i < outcomes_files_size; i++) {

                    for (int t = 0; t < input_params.tier_depth; t++) {

                        fprintf(
                            outcomes_file[i], 
                            "%d,%-.8E,%d,%-.8E,%s,%s,%-.8E\n", 
                            t, input_params.fp_prob[t],
                            t, input_params.alpha[t],
                            PENALTY_TYPE_STR[input_params.fp_resolution_tech],  
                            (*leaf_itr)->get_outcome().c_str(), 
                            (*leaf_itr)->get_prob_val());
                    }

                    // fprintf(
                    //     input_params.alpha_outcomes_file[i], 
                    //     "%-.8E,%s,%-.8E\n", 
                    //     input_params.alpha[i], (*leaf_itr)->get_outcome().c_str(), (*leaf_itr)->get_prob_val());
                }
            }
        }

        ++leaf_itr;
    }

    if ((modes & MODE_VERBOSE)) {

        printf("\n[CHECKSUM : %-.8E]\n", checksum);
        printf("[AVG_LATENCY : %-.8E]\n", avg_latency);
        printf("[CACHE_LATENCY : %-.8E]\n", latency_to_content(input_params.content_sources, input_params.latencies));

        double origin_latency = latency_through_tier(input_params.origin_tier, input_params.latencies);
        printf("[ORIGIN_LATENCY : %-.8E]\n", origin_latency);
    }

    // erase the whole tree...
    erase_decision_tree(decision_tree);

    return 0;
}

