#include <math.h>
#include <cfloat>
#include <algorithm>

#include "tree/tree.hh"
#include "node.h"
#include "graph.h"
#include "rid-analytics.h"

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
    for (int t = 1; t < tier; t++) {
        penalty_candidate += 2.0 * latencies[t];
    }

    // the 'peer' link latency part
    // FIXME: this isn't always correct but let it stay like that for now
    penalty_candidate = latencies[tier];

    return penalty_candidate + base_latency;
}

double latency_to_content(int * content_sources, double * latencies) {

    double latency_to_content = 0.0 * latencies[0];

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

void calc_penalties(struct RIDAnalytics::rid_analytics_inputs & input_params) {

    input_params.penalties = (double *) calloc(input_params.tier_depth, sizeof(double));

    // first tier at which content is available
    int origin_tier = 0;
    for (int i = 0; i < input_params.tier_depth; i++) {

        if (input_params.content_sources[i] > 0) {
            
            origin_tier = i;
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
            input_params.penalties[m] += latency_through_tier(origin_tier, input_params.latencies);
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

                for (int i = m; i >= m; i--) {
                    
                    // probabilities always < 1.0 please: if there are more 
                    // content sources than domains (and assuming a uniform 
                    // distribution of sources per domains), we truncate 
                    // probabilities here...
                    prob_next = min(((double) input_params.content_sources[m]) / ((double) (input_params.domains[i])), 1.0);
                    amortized_penalty += prob_prev * (1.0 - prob_next) * latency_through_tier(m, input_params.latencies);

                    prob_prev = prob_prev * prob_next;
                }

                // special operation for the minimum latency case...
                // FIXME: not sure if totally correct
                amortized_penalty += prob_prev * 2.0 * input_params.latencies[0];
            }

            input_params.penalties[m] = amortized_penalty;
        }
    }
}

void print_params(struct RIDAnalytics::rid_analytics_inputs input_params) {

    printf("\n*** PARAMETERS ***\n");

    // initial row w/ tier numbers
    printf("\n\n[TIER]    : |");

    for (int i = 0; i < input_params.tier_depth; i++)
        printf(" %-8d|", i + 1);

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
    std::string data_dir) {

    // latency values per tier, in multiples of some time unit.
    double latencies[MAX_TIER_DEPTH];
    // FIXME: just initialize everything to 1.0
    for (int tier = 0; tier < input_params.tier_depth; tier++)
        latencies[tier] = 1.0;

    // point the latencies array in the input params struct to the array 
    // just created
    input_params.latencies = latencies;

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
    Node * root_node = new Node(0, 0, 0, 1.0, 0.0, Node::ORI_NODE);

    root = decision_tree.begin();
    root = decision_tree.insert(root, root_node);

    // we fill the tree tier by tier. we add children to each non-leaf node 
    // of the previous tier. to iterate over nodes of a given tier, 
    // we use fixed_depth_iterators of tree.hh.
    tree<Node *>::fixed_depth_iterator depth_itr;
    tree<Node *>::fixed_depth_iterator end_depth_itr;

    Node * mhs_node = NULL;
    Node * mhd_node = NULL;
    Node * sfp_node = NULL;
    Node * tpo_node = NULL;
    Node * def_node = NULL;

    int curr_node_tier = 0;
    int curr_node_max_tier = 0;
    double prev_node_prob = 0.0;
    double prev_node_latency = 0.0;
    Node::Type prev_node_type = Node::UNKNOWN;

    double path_latency = 0.0;
    double penalty_latency = 0.0;

    double prob_fp_here = 0.0;
    double p = 0.0;

    // start the .dot file for rendering the probability tree
    Graph * graphviz_graph = NULL;
    if ((modes & MODE_SAVEGRAPH))
        graphviz_graph = new Graph(std::string(data_dir + "/graphviz").c_str());

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

            // save the information about the node to which we will append the 
            // C, I and N children
            prev_node_prob = (*depth_itr)->get_prob_val();
            prev_node_latency = (*depth_itr)->get_latency_val();
            prev_node_type = (*depth_itr)->get_type();

            // create the children Node objects, specifying the diff types
            mhs_node = new Node(curr_node_max_tier, curr_node_tier, Node::MHS_NODE);
            mhd_node = new Node(curr_node_max_tier, curr_node_tier, Node::MHD_NODE);
            tpo_node = new Node(curr_node_max_tier, curr_node_tier, Node::TPO_NODE);
            sfp_node = new Node(curr_node_max_tier, curr_node_tier, Node::SFP_NODE);
            def_node = new Node(curr_node_max_tier, curr_node_tier, Node::DEF_NODE);

            if (prev_node_type == Node::ORI_NODE || prev_node_type == Node::DEF_NODE) {

                // Q3) which TP reality are we in ? no TPs? single TP? 
                // multiple TPs? 
                if (input_params.content_sources[curr_node_tier] == 0) {         // TP reality = no TPs

                    // *** probabilities ***
                    // FIXME: p is completely incorrect
                    p = 0.001 * input_params.fp_prob[curr_node_tier];
                    prob_fp_here = (input_params.fp_prob[curr_node_tier] - p);
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(p * prev_node_prob);
                    def_node->set_prob_val((1.0 - input_params.fp_prob[curr_node_tier]) * prev_node_prob);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    mhd_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    tpo_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    sfp_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    def_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(curr_node_tier);
                    mhd_node->set_next_tier(curr_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(curr_node_tier);

                    // default route : go up a tier!
                    def_node->set_next_tier(curr_node_tier + 1);
                    def_node->set_max_tier(curr_node_tier + 1);

                    if ((modes & MODE_SAVEGRAPH)) {

                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);
                        graphviz_graph->add_node((*depth_itr), depth, sfp_node, depth + 1, breadth, p * prev_node_prob);
                        graphviz_graph->add_node((*depth_itr), depth, def_node, depth + 1, breadth, (1.0 - input_params.fp_prob[curr_node_tier]));

                        Node * nodes[] = {mhs_node, mhd_node, sfp_node, def_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] == 1) {   // TP reality = 1 TP

                    // *** probabilities ***
                    prob_fp_here = input_params.fp_prob[curr_node_tier];
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val((1.0 - input_params.fp_prob[curr_node_tier]) * prev_node_prob);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    mhd_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    tpo_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    sfp_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    def_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(curr_node_tier);
                    mhd_node->set_next_tier(curr_node_tier);
                    tpo_node->set_next_tier(curr_node_tier);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        // add nodes to the .dot file for graph rendering
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);
                        graphviz_graph->add_node((*depth_itr), depth, tpo_node, depth + 1, breadth, (1.0 - input_params.fp_prob[curr_node_tier]));

                        Node * nodes[] = {mhs_node, mhd_node, tpo_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] > 1) {    // TP reality = n TPs

                    // *** probabilities ***
                    mhs_node->set_prob_val(input_params.alpha[curr_node_tier] * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    mhd_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    tpo_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    sfp_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);
                    def_node->set_latency_val(prev_node_latency + latencies[curr_node_tier]);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(curr_node_tier);
                    mhd_node->set_next_tier(curr_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, input_params.alpha[curr_node_tier]);
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]));

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::SFP_NODE) {

                // Q1) if the previous decision was SFP_NODE, we will now go 
                // down a tier
                int next_node_tier = curr_node_tier - 1;

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_tier == END_OF_PATH) {

                    prob_fp_here = input_params.fp_prob[curr_node_tier];
                    
                    // penalty due to relaying
                    path_latency = 0.0;
                    penalty_latency = input_params.penalties[curr_node_max_tier];

                    mhd_node->set_outcome(std::string(OUTCOME_DROPPED));
                    sfp_node->set_outcome(std::string(OUTCOME_IDEST_CSERVER));

                } else {

                    prob_fp_here = input_params.fp_prob[next_node_tier];

                    // the next component of path latency will be that of 
                    // the tier below: notice that we're essentially forwarding 
                    // between hops within the same tier (in this case, 
                    // next_node_tier)
                    path_latency = latencies[next_node_tier];
                    penalty_latency = 0.0;
                } 

                // Q3) which TP reality are we in ? no TPs? single TP? 
                // multiple TPs? 
                if (input_params.content_sources[curr_node_tier] == 0) {         // TP reality = no TPs

                    // *** probabilities ***
                    // FIXME: p is completely incorrect
                    mhs_node->set_prob_val((prob_fp_here) * input_params.alpha[curr_node_tier] * prev_node_prob);
                    mhd_node->set_prob_val((prob_fp_here) * (1.0 - input_params.alpha[curr_node_tier]) * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(1.0 - prob_fp_here);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(next_node_tier);
                    def_node->set_next_tier(END_OF_PATH);

                    mhs_node->set_outcome(std::string(OUTCOME_IDEST_CSERVER));

                    if ((modes & MODE_SAVEGRAPH)) {

                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (prob_fp_here) * input_params.alpha[curr_node_tier]);
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (prob_fp_here) * (1.0 - input_params.alpha[curr_node_tier]));
                        graphviz_graph->add_node((*depth_itr), depth, sfp_node, depth + 1, breadth, 1.0 - prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node, sfp_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] == 1) {   // TP reality = 1 TP

                    // *** probabilities ***
                    prob_fp_here = 1.0;
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    mhs_node->set_outcome(std::string(OUTCOME_CCACHE));

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        // add nodes to the .dot file for graph rendering
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] > 1) {    // TP reality = n TPs

                    // *** probabilities ***
                    prob_fp_here = 1.0;
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    mhs_node->set_outcome(std::string(OUTCOME_CCACHE));

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        // add nodes to the .dot file for graph rendering
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::TPO_NODE) {

                // Q1) if the previous decision was TPO_NODE, we will now go 
                // down a tier
                int next_node_tier = curr_node_tier - 1;

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_tier == END_OF_PATH) {

                    // penalty due to relaying
                    path_latency = 0.0;
                    penalty_latency = input_params.penalties[curr_node_max_tier];

                    mhs_node->set_outcome(std::string(OUTCOME_CCACHE));
                    mhd_node->set_outcome(std::string(OUTCOME_DROPPED));
                    tpo_node->set_outcome(std::string(OUTCOME_CCACHE));

                } else {

                    // the next component of path latency will be that of 
                    // the tier below: notice that we're essentially forwarding 
                    // between hops within the same tier (in this case, 
                    // next_node_tier)
                    path_latency = latencies[next_node_tier];
                    penalty_latency = 0.0;
                } 

                // Q3) which TP reality are we in ? no TPs? single TP? 
                // multiple TPs? 
                if (input_params.content_sources[curr_node_tier] == 0) {         // TP reality = no TPs

                    // *** probabilities ***
                    mhs_node->set_prob_val(0.0);
                    mhd_node->set_prob_val(0.0);
                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(END_OF_PATH);
                    mhd_node->set_next_tier(END_OF_PATH);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                } else if (input_params.content_sources[curr_node_tier] == 1) {   // TP reality = 1 TP

                    // *** probabilities ***
                    prob_fp_here = input_params.fp_prob[curr_node_tier];
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val((1.0 - prob_fp_here) * prev_node_prob);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(next_node_tier);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, 1.0 - prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node, tpo_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] > 1) {    // TP reality = n TPs

                    // *** probabilities ***
                    prob_fp_here = 1.0;
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {
                        
                        // add nodes to the .dot file for graph rendering
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::MHS_NODE) {

                // Q1) if the previous decision was MHS_NODE, we will now go 
                // down a tier
                int next_node_tier = curr_node_tier - 1;

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_tier == END_OF_PATH) {

                    // penalty due to relaying
                    path_latency = 0.0;
                    penalty_latency = input_params.penalties[curr_node_max_tier];

                    mhs_node->set_outcome(std::string(OUTCOME_CCACHE));
                    mhd_node->set_outcome(std::string(OUTCOME_DROPPED));

                } else {

                    // the next component of path latency will be that of 
                    // the tier below: notice that we're essentially forwarding 
                    // between hops within the same tier (in this case, 
                    // next_node_tier)
                    path_latency = latencies[next_node_tier];
                    penalty_latency = 0.0;
                } 

                // Q3) which TP reality are we in ? no TPs? single TP? 
                // multiple TPs? 
                if (input_params.content_sources[curr_node_tier] == 0) {         // TP reality = no TPs

                    // *** probabilities ***
                    // FIXME: p is completely incorrect
                    p = 0.001 * input_params.fp_prob[curr_node_tier];
                    prob_fp_here = (input_params.fp_prob[curr_node_tier] - p);
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {

                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] == 1) {   // TP reality = 1 TP

                    // *** probabilities ***
                    prob_fp_here = (input_params.fp_prob[curr_node_tier]);
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {

                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else if (input_params.content_sources[curr_node_tier] > 1) {    // TP reality = n TPs

                    // *** probabilities ***
                    prob_fp_here = 1.0;
                    mhs_node->set_prob_val((input_params.alpha[curr_node_tier] * prob_fp_here) * prev_node_prob);
                    mhd_node->set_prob_val((1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here * prev_node_prob);

                    tpo_node->set_prob_val(0.0);
                    sfp_node->set_prob_val(0.0);
                    def_node->set_prob_val(0.0);

                    // *** latencies *** 
                    mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                    def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                    // *** next tier *** 
                    mhs_node->set_next_tier(next_node_tier);
                    mhd_node->set_next_tier(next_node_tier);
                    tpo_node->set_next_tier(END_OF_PATH);
                    sfp_node->set_next_tier(END_OF_PATH);
                    def_node->set_next_tier(END_OF_PATH);

                    if ((modes & MODE_SAVEGRAPH)) {

                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));
                        graphviz_graph->add_node((*depth_itr), depth, mhd_node, depth + 1, breadth, (1.0 - input_params.alpha[curr_node_tier]) * prob_fp_here);

                        Node * nodes[] = {mhs_node, mhd_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::MHD_NODE) {

                // Q1) if the previous decision was MHD_NODE, we will ALWAYS 
                // end the path

                // penalty due to relaying
                path_latency = 0.0;
                penalty_latency = input_params.penalties[curr_node_max_tier];

                mhd_node->set_outcome(std::string(OUTCOME_DROPPED));

                mhs_node->set_prob_val(0.0);
                mhd_node->set_prob_val(prev_node_prob * 1.0);
                tpo_node->set_prob_val(0.0);
                sfp_node->set_prob_val(0.0);
                def_node->set_prob_val(0.0);

                // *** latencies *** 
                mhs_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                mhd_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                tpo_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                sfp_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                def_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);                    

                // *** next tier *** 
                mhs_node->set_next_tier(END_OF_PATH);
                mhd_node->set_next_tier(END_OF_PATH);
                tpo_node->set_next_tier(END_OF_PATH);
                sfp_node->set_next_tier(END_OF_PATH);
                def_node->set_next_tier(END_OF_PATH);

                if ((modes & MODE_SAVEGRAPH)) {

                    // add nodes to the .dot file for graph rendering
                    graphviz_graph->add_node((*depth_itr), depth, mhs_node, depth + 1, breadth, (input_params.alpha[curr_node_tier] * prob_fp_here));

                    Node * nodes[] = {mhd_node};
                    graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                }
            }

            // finally append the 3 children
            decision_tree.append_child(depth_itr, mhs_node);
            decision_tree.append_child(depth_itr, mhd_node);
            decision_tree.append_child(depth_itr, sfp_node);
            decision_tree.append_child(depth_itr, tpo_node);
            decision_tree.append_child(depth_itr, def_node);

            breadth++;

            // go to the next node in the previous tier and repeat the process
            ++depth_itr;
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
                // [FP PROB AT TIER i], [LOOKUP OUTCOME TYPE], [OUTCOME TYPE PROBABILITY OF OCCURRENCE]
                for (int i = 0; i < input_params.tier_depth; i++) {

                    fprintf(
                        outcomes_file[i], 
                        "%-.8E,%s,%-.8E\n", 
                        input_params.fp_prob[i], (*leaf_itr)->get_outcome().c_str(), (*leaf_itr)->get_prob_val());

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

        double origin_latency = latency_through_tier(input_params.origin_tier, latencies);
        printf("[ORIGIN_LATENCY : %-.8E]\n", origin_latency);
    }

    // erase the whole tree...
    erase_decision_tree(decision_tree);

    return 0;
}
