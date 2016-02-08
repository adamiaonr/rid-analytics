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

        if ((*leaf_itr)->get_next_level() != END_OF_PATH) {

            return 0;
        }

        ++leaf_itr;
    }

    return 1;
}

double latency_to_level(int level, double * latencies) {

    // go up to level, then down
    double penalty_candidate = 0.0 * latencies[0] + latencies[level];

    for (int l = 0; l < level; l++) {
        penalty_candidate += 2.0 * latencies[l];
    }

    return penalty_candidate;
}

double latency_to_content(double * latencies, int cache_level) {

    double latency_to_content = 0.0 * latencies[0];

    for (int l = 0; l < MAX_LEVEL_DEPTH; l++) {

         if (l == cache_level) {
            
            latency_to_content += latencies[l];

            // 'that's it, the rebels are there!'
            break;
         }

         latency_to_content += 2 * latencies[l];
    }

    return latency_to_content;
}

int RIDAnalytics::run_model(
    double & avg_latency,
    double * fp_prob, 
    double * o_optimistic, 
    int level_depth,
    int cache_level, 
    int origin_level,
    int fp_resolution_tech,
    bool verbose,
    bool save_cdf,
    bool save_graph,
    bool save_outcomes,
    FILE ** fp_prob_outcomes_file,
    FILE ** o_optimistic_outcomes_file,
    std::string data_dir) {

    // latency values per level, in multiples of some time unit.
    double latencies[MAX_LEVEL_DEPTH];
    // FIXME: just initialize everything to 1.0
    for (int lvl = 0; lvl < level_depth; lvl++)
        latencies[lvl] = 1.0;

    // the penalty array. our goal is to calculate the relaying 
    // penalty - PENALTY(m) - for a request, given the following conditions: 
    //
    //  -# the request has crossed a level as high as level m
    //  -# we know the distribution of sources per level for the scenario 
    //     being evaluated
    //  -# we know the avg. latencies of forwarding between routers at level l 
    //
    // this can be pre-calculated 'a priori'. how? see below...
    //
    double penalties[MAX_LEVEL_DEPTH];

    double penalty_candidate = DBL_MAX;
//    double penalty_min = DBL_MAX;
//    double penalty_max = 0.0;

    if (fp_resolution_tech == PENALTY_TYPE_FEEDBACK) {

        // this case is straight forward: if we're at a destination after 
        // passing through level m, the penalty will be composed by the 
        // following parts:
        //
        //  1) latency of going back to the source of the request, i.e. 
        //      2 * L[0] + 2 * L[1] + ... + 2 * L[m - 1] + L[m] 
        //  2) latency of going up to the first visible origin server
        //      2 * L[0] + 2 * L[1] + ... + 2 * L[v - 1] + L[v]
        //

        // calculate a penalty for each level max. m
        for (int m = 0; m < level_depth; m++) {

            // 1) going back from m to content source
            penalties[m] = latency_to_level(m, latencies);

            // 2) we got to the source of the request, now moving back up, till 
            // a content source is visible
            //penalties[m] += latency_to_content(latencies, content_source_dist);
            penalties[m] += latency_to_level(origin_level, latencies);
        }

    } else if (fp_resolution_tech == PENALTY_TYPE_FALLBACK) {

        for (int m = 0; m < level_depth; m++) {

            // FIXME: legacy is stupid...
            int l = origin_level;

            // say that we can go up to level m and get to a visible 
            // content source. how long would that take?
            if (l > m) {

                // easy: go up to level l, then down
                penalty_candidate = latency_to_level(l, latencies);

            } else if (l == m) {

                // say there are content sources that become visible at 
                // our level. we may even be 3 hops away from them! what's 
                // the probability of that happening?
                // at level l we have content_source_dist[l] sources that 
                // become visible. these can be distributed over 
                // breadth[l] * breadth[l - 1] * ... * breadth[0] different 
                // slots.
                double slots = 1.0;
                double level_breadth = 4.0;

                for (int a = l; a >= l; a--) {
                    slots *= slots * (double) level_breadth;
                }

                // therefore, P('3 hops') = cs[l] / slots
                double prob_3_hops = min((double) ((double) 1 / slots), 1.0);

                // besides the '3 hops', we can have the penalty of going 
                // to another domain of level l.
                double the_other = latency_to_level(l, latencies);

                // the penalty will then be the weighted average of these
                // two latencies. this is the minimum value we'll ever get.
                penalty_candidate = prob_3_hops * (1.0 * latencies[0]) +
                    (1.0 - prob_3_hops) * the_other;

            } else if (l < m) {

                // say we have to go down a level to get visible sources. 
                // we first go up to m, then go down to 0.
                penalty_candidate = latency_to_level(l, latencies);
            }

            // the FALLBACK penalty for level max. m becomes the minimum from 
            // all values gathered above.
            penalties[m] = penalty_candidate;
        }
    }

    // print all parameters (if verbose mode is set)
    if (verbose) {

        printf("\n*** PARAMETERS ***\n");

        // initial row & midrule
        printf("\n\n[LEVEL]   : |");

        for (int i = 0; i < level_depth; i++)
            printf(" %-8d|", i + 1);

        printf("\n-------------");

        for (int i = 0; i < level_depth; i++)
            printf("----------");

        // penalty row per level
        printf("\n[PENALTY] : |");

        for (int i = 0; i < level_depth; i++)
            printf(" %-.6f|", penalties[i]);

        printf("\n[FP-PROB] : |");

        for (int i = 0; i < level_depth; i++)
            printf(" %-.2E|", fp_prob[i]);

        printf("\n[O-OPT]   : |");

        for (int i = 0; i < level_depth; i++)
            printf(" %-.2E|", o_optimistic[i]);

        printf("\n\n");
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
    Node * root_node = new Node(0, 0, 0, 1.0, 0.0, Node::O_NODE);
    root = decision_tree.begin();
    root = decision_tree.insert(root, root_node);

    // we fill the tree level by level. we add children to each non-leaf node 
    // of the previous level. to iterate over nodes of a given level, 
    // we use fixed_depth_iterators of tree.hh.
    tree<Node *>::fixed_depth_iterator depth_itr;
    tree<Node *>::fixed_depth_iterator end_depth_itr;

    // this is an iterative process, each iteration corresponding to a 
    // forwarding decision at a router. our model considers 3 
    // possible forwarding decisions:
    //  -# (I)ncorrect : request is forwarded over a 'wrong' interface, 
    //     i.e. picking a FP entry which doesn't point to the same 
    //     interface as the TP entry (if any)
    //  -# (C)orrect : the request is sent over the same interface as the 
    //     TP entry. this can happen when the TP is picked OR when a FP 
    //     pointing to the same interface as the TP is picked.
    //  -# (N)eutral : no matches in the forwarding table, a TN occurs
    //
    Node * i_node = NULL;
    Node * c_node = NULL;
    Node * n_node = NULL;

    int curr_node_level = 0;
    int curr_node_max_level = 0;
    double prev_node_prob = 0.0;
    double prev_node_latency = 0.0;
    Node::Type prev_node_type = Node::UNKNOWN;

    double path_latency = 0.0;
    double penalty_latency = 0.0;

    // start the .dot file for rendering the probability tree
    Graph * graphviz_graph = NULL;
    if (save_graph)
        graphviz_graph = new Graph(std::string(data_dir + "/graphviz").c_str());

    // depth and breadth are sort of (x, y) coordinates for nodes in a tree. 
    // these will be used to generate node names for the .dot file. the origin 
    // of the tree has coordinates (0, 0).
    int depth = 0;
    int breadth = 0;

    // note the limit condition : on the longest possible path (i.e. going 
    // through the max. level) we make 2 * level_depth decisions
    for (depth = 0; depth < ((2 * level_depth)); depth++) {

        // update the iterators
        depth_itr = decision_tree.begin_fixed(root, depth);

        breadth = 0;

        // cycle through the previous decisions, append children
        while (decision_tree.is_valid(depth_itr)) {

            curr_node_level = (*depth_itr)->get_next_level();
            curr_node_max_level = (*depth_itr)->get_max_level();

            // if this is a 'dead end', don't append children to it
            if (curr_node_level == END_OF_PATH) {
                
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
            i_node = new Node(curr_node_max_level, curr_node_level, Node::I_NODE);
            c_node = new Node(curr_node_max_level, curr_node_level, Node::C_NODE);
            n_node = new Node(curr_node_max_level, curr_node_level, Node::N_NODE);

            // what was the previous decision type? this will influence the 
            // next decision: 
            // 
            //  PREV   | NEXT
            // --------------------
            //  N or O | I, C or N
            //  I      | I or N
            //  C      | I or C       
            //
            // for each NEXT decision, we ask/answer 3 questions:
            //  Q1) at what level will the NEXT decision leave the request?
            //  Q2) what latency does the NEXT decision add?
            //  Q3) what is the probability of making the NEXT decision?
            //
            // if a decision type does not show up in the 'NEXT' column, it is 
            // given a prob of 0.0 and it's next level is set to END_OF_PATH.  
            //
            if (prev_node_type == Node::N_NODE || prev_node_type == Node::O_NODE) {

                // Q1) if NEXT is N : +1 level (and update the max level)
                n_node->set_next_level(curr_node_level + 1);
                n_node->set_max_level(curr_node_level + 1);
                // Q1) if NEXT is I or C : next level stays the same
                i_node->set_next_level(curr_node_level);
                c_node->set_next_level(curr_node_level);

                // Q2) for all decisions, we add +1 latency of the curr_node_level
                i_node->set_latency_val(prev_node_latency + latencies[curr_node_level]);
                c_node->set_latency_val(prev_node_latency + latencies[curr_node_level]);
                n_node->set_latency_val(prev_node_latency + latencies[curr_node_level]);

                // Q3) the probability depends on the content sources which 
                // are 'visible' at the current level
                if (curr_node_level == cache_level) {

                    // if there are 'visible' sources, then there's no way we 
                    // can get an N decision. mark it as END_OF_PATH.
                    n_node->set_prob_val(prev_node_prob * 0.0);
                    n_node->set_next_level(END_OF_PATH);

                    // we may then follow a C or I decision. as mentioned above, 
                    // the probability of an incorrect decision is only due to 
                    // FPs, multiplied by the 'pessimistic' factor.
                    double prob_fpi = (1.0 - o_optimistic[curr_node_level]) * fp_prob[curr_node_level]; 
                    
                    // correct decisions can happen due to 'optimistic' FPs and 
                    // TPs. note that in this case (in which TNs are 
                    // impossible) 1 - P(FP) corresponds to the probability of 
                    // having a router follow the TP. therefore, the 
                    // probability of a correct decision is 1.0 - prob_fpi.
                    i_node->set_prob_val(prev_node_prob * prob_fpi);
                    c_node->set_prob_val(prev_node_prob * (1.0 - prob_fpi));

                    if (save_graph) {
                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, prob_fpi);
                        graphviz_graph->add_node((*depth_itr), depth, c_node, depth + 1, breadth, 1.0 - prob_fpi);

                        Node * nodes[] = {i_node, c_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }

                } else {

                    // if there *NO* 'visible' sources, then there's no chance 
                    // of making a C decision: all FPs will be (I)ncorrect. 
                    // mark it as END_OF_PATH.
                    c_node->set_prob_val(prev_node_prob * 0.0);
                    c_node->set_next_level(END_OF_PATH);

                    // we may then follow an I or N decision, which basically 
                    // comes down to P(FP) or (1 - P(FP)).
                    // *********************************************************
                    // FIXME: should we add the o_optimistic value here?
                    // *********************************************************
                    n_node->set_prob_val(prev_node_prob * (1.0 - fp_prob[curr_node_level]));
                    i_node->set_prob_val(prev_node_prob * fp_prob[curr_node_level]);

                    if (save_graph) {

                        // add nodes to the .dot file for graph rendering
                        graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, fp_prob[curr_node_level]);
                        graphviz_graph->add_node((*depth_itr), depth, n_node, depth + 1, breadth, (1.0 - fp_prob[curr_node_level]));

                        Node * nodes[] = {i_node, n_node};
                        graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                    }
                }

            } else if (prev_node_type == Node::I_NODE) {

                // if the previous decision was I, we can either make another I 
                // decision or come to a dead e(N)d (N will require relaying as 
                // a solution).

                // making a C decision is impossible
                c_node->set_prob_val(prev_node_prob * 0.0);

                // mark the C and N nodes as END_OF_PATH.
                c_node->set_next_level(END_OF_PATH);
                n_node->set_next_level(END_OF_PATH);
                n_node->set_outcome(std::string(OUTCOME_DROPPED));

                // Q1) if the previous decision was I, we will now go down 
                // a level
                int next_node_level = curr_node_level - 1;
                double prob_incorrect = 0.0;

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_level == END_OF_PATH) {

                    // penalty due to relaying
                    path_latency = 0.0;
                    penalty_latency = penalties[curr_node_max_level];

                    i_node->set_outcome(std::string(OUTCOME_IDEST_CSERVER));

                    // we did not change levels, and so fwd tables will remain 
                    // the same. we'll make the same decision again.
                    prob_incorrect = 1.0;

                } else {

                    // the next component of path latency will be that of 
                    // the level below: notice that we're essentially forwarding 
                    // between hops within the same level (in this case, 
                    // next_node_level)
                    path_latency = latencies[next_node_level];
                    penalty_latency = 0.0;

                    prob_incorrect = fp_prob[next_node_level];
                }   

                i_node->set_next_level(next_node_level);

                // Q2) the P(FP) we use is that of the next level
                n_node->set_prob_val(prev_node_prob * (1.0 - prob_incorrect));
                i_node->set_prob_val(prev_node_prob * prob_incorrect);

                // Q3) regarding latencies:
                //  I node : use the values calculated above
                i_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                //  N node : the policy is to relay when a I > N transition 
                //           happens, so simply add this level's penalty
                n_node->set_latency_val(prev_node_latency + penalties[curr_node_max_level]);

                if (save_graph) {

                    // add nodes to the .dot file for graph rendering
                    graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, prob_incorrect);
                    graphviz_graph->add_node((*depth_itr), depth, n_node, depth + 1, breadth, (1.0 - prob_incorrect));

                    Node * nodes[] = {i_node, n_node};
                    graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                }

            } else {

                // if the previous decision was C, we can either make another 
                // C or an I decision. mark the N node as END_OF_PATH.
                n_node->set_prob_val(prev_node_prob * 0.0);
                n_node->set_next_level(END_OF_PATH);

                // Q1) if the previous decision was C, we will also go down 
                // a level
                int next_node_level = curr_node_level - 1;
                double prob_incorrect = 0.0;

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_level == END_OF_PATH) {

                    // penalty due to relaying
                    path_latency = 0.0;
                    penalty_latency = penalties[curr_node_max_level];

                    c_node->set_outcome(std::string(OUTCOME_CCACHE));
                    i_node->set_outcome(std::string(OUTCOME_IDEST_CSERVER));

                    // notice this is the probability of an incorrect decision
                    prob_incorrect = (1.0 - o_optimistic[curr_node_level]) * fp_prob[curr_node_level];

                } else {

                    path_latency = latencies[next_node_level];
                    penalty_latency = 0.0;

                    prob_incorrect = (1.0 - o_optimistic[next_node_level]) * fp_prob[next_node_level];
                }

                i_node->set_next_level(next_node_level);
                c_node->set_next_level(next_node_level);

                // Q2) probabilities 
                i_node->set_prob_val(prev_node_prob * prob_incorrect);
                c_node->set_prob_val(prev_node_prob * (1.0 - prob_incorrect));

                // Q3) latencies
                i_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                c_node->set_latency_val(prev_node_latency + path_latency);

                if (save_graph) {

                    // add nodes to the .dot file for graph rendering
                    graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, prob_incorrect);
                    graphviz_graph->add_node((*depth_itr), depth, c_node, depth + 1, breadth, (1.0 - prob_incorrect));

                    Node * nodes[] = {i_node, c_node};
                    graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
                }
            }

            // finally append the 3 children
            decision_tree.append_child(depth_itr, i_node);
            decision_tree.append_child(depth_itr, c_node);
            decision_tree.append_child(depth_itr, n_node);

            breadth++;

            // go to the next node in the previous level and repeat the process
            ++depth_itr;
        }

        if (is_all_eop(decision_tree) > 0)
            break;
    }

    // wrap the .dot file up
    if (save_graph)
        graphviz_graph->terminate();

    // the decision tree is finished. let's iterate through the leafs to learn 
    // a bit about latencies for this scenario.
    tree<Node *>::leaf_iterator leaf_itr = decision_tree.begin_leaf();
    tree<Node *>::leaf_iterator end_leaf_itr = decision_tree.end_leaf();

    // print some results for quick checking
    if (verbose)
        printf("\n*** FINAL RESULTS ***\n\n");

    double cache_latency = DBL_MAX;
    double checksum = 0.0;
    
    char * node_str = NULL;

    while (leaf_itr != end_leaf_itr) {

        if ((*leaf_itr)->get_prob_val() > 0.0) {

            if (verbose) {

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

            if (save_outcomes) {

                // if requested, save outcome information in given files
                for (int i = 0; i < level_depth; i++) {

                    fprintf(
                        fp_prob_outcomes_file[i], 
                        "%-.8E,%s,%-.8E\n", 
                        fp_prob[i], (*leaf_itr)->get_outcome().c_str(), (*leaf_itr)->get_prob_val());

                    fprintf(
                        o_optimistic_outcomes_file[i], 
                        "%-.8E,%s,%-.8E\n", 
                        o_optimistic[i], (*leaf_itr)->get_outcome().c_str(), (*leaf_itr)->get_prob_val());
                }
            }
        }

        ++leaf_itr;
    }

    if (verbose) {

        printf("\n[CHECKSUM : %-.8E]\n", checksum);
        printf("[AVG_LATENCY : %-.8E]\n", avg_latency);
        printf("[CACHE_LATENCY : %-.8E]\n", latency_to_content(latencies, cache_level));

        double origin_latency = latency_to_level(origin_level, latencies);
        printf("[ORIGIN_LATENCY : %-.8E]\n", origin_latency);
    }

    // erase the whole tree...
    erase_decision_tree(decision_tree);

    return 0;
}
