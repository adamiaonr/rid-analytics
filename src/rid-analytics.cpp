#include <math.h>
#include <cfloat>
#include <algorithm>

#include "argvparser.h"
#include "tree/tree.hh"

#include "node.h"
#include "dataparser.h"
#include "graph.h"

#define MAX_FILENAME_SIZE   64
#define MAX_REQUEST_SIZE    30
#define MAX_LEVEL_DEPTH     30

#define DEFAULT_O_OPTIMISTIC (double) (1.0 / 10.0)
#define END_OF_PATH (int) -1

#define OPTION_SCN_FILE         (char *) "scn-file"
#define OPTION_SCN_DIR          (char *) "scn-dir"
#define OPTION_DAT_DIR          (char *) "dat-dir"
#define OPTION_GRAPHVIZ_FILE    (char *) "graphviz-file"

#define CACHING             (char *) "caching"
#define LEVEL_DEPTH         (char *) "level_depth"
#define LEVEL_BREADTH       (char *) "level_breadth"
#define CONTENT_SOURCE_DIST (char *) "content_source_dist"
#define COMPLEMENT_DIST     (char *) "complement_dist"
#define LATENCIES           (char *) "latencies"
#define REQUEST_SIZE        (char *) "request_size"

#define BF_SIZE 160

using namespace std;
using namespace CommandLineProcessing;

enum DecisionMode {
    DEFAULT = 0x00
};

ArgvParser * create_argv_parser() {

    ArgvParser * cmds = new ArgvParser();

    cmds->setIntroductoryDescription("\n\nrid-analytics v0.1\n\ntool to run simple analytical "\
        "evaluations on networks which use RIDs (i.e. Bloom Filters) for packet "\ 
        "forwarding. it computes probability distributions for latencies and "\
        "avg. latencies/ relaying penalties given a set of scenario characteristics (e.g. use of "\
        "caching, distribution of content sources, etc.).\nby adamiaonr@cmu.edu");

    cmds->setHelpOption("h", "help", "help page.");

    cmds->defineOption(
            OPTION_SCN_FILE,
            "path to .scn file which contains info about scenario to evaluate",
            ArgvParser::OptionRequiresValue);

    // cmds->defineOption(
    //         OPTION_SCN_DIR,
    //         "path to a directory with multiple .scn files to evaluate",
    //         ArgvParser::OptionRequiresValue);

    cmds->defineOption(
            OPTION_DAT_DIR,
            "path to a directory with .dat files with |F\\R| distributions",
            ArgvParser::OptionRequiresValue);

        cmds->defineOption(
            OPTION_GRAPHVIZ_FILE,
            "name of the graphviz file to output (will generate a <filename>.dot file)",
            ArgvParser::OptionRequiresValue);

    return cmds;
}

double fp_rate(double m, double n, double k, double c) {

    return pow((1.0 - exp(-((n / m) * k))), k * c);
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

int main (int argc, char **argv) {

    // ************************************************************************
    // 1) initialization: arguments
    // ************************************************************************

    ArgvParser * cmds = create_argv_parser();

    // paths for .scn file and/or directory of .scn files 
    char * scn_file;
    //char * scn_dir;
    char * graphviz_file;

    // path for .dat file dir
    char * dat_dir;

    // trigger argument parsing with parse()
    int result = cmds->parse(argc, argv);

    if (result != ArgvParser::NoParserError) {

        // something went wrong. show help option.
        fprintf(stderr, "%s\n", cmds->parseErrorDescription(result).c_str());

        if (result != ArgvParser::ParserHelpRequested) {
            fprintf(stderr, "use option -h for help.\n");
        }

        return -1;

    } else {

        if (cmds->foundOption(OPTION_SCN_FILE)) {

            scn_file = (char *) cmds->optionValue(OPTION_SCN_FILE).c_str();

        } else {

            fprintf(stderr, "no scenario file path specified. use "\
                "option -h for help.\n");

            return -1;
        }

        if (cmds->foundOption(OPTION_DAT_DIR)) {

            dat_dir = (char *) cmds->optionValue(OPTION_DAT_DIR).c_str();

        } else {

            fprintf(stderr, "no .dat dir path specified. use "\
                "option -h for help.\n");

            return -1;
        }

        if (cmds->foundOption(OPTION_GRAPHVIZ_FILE)) {

            graphviz_file = (char *) cmds->optionValue(OPTION_GRAPHVIZ_FILE).c_str();

        } else {

            fprintf(stderr, "no graphviz file name specified. use "\
                "option -h for help.\n");

            return -1;
        }

/*        if (cmds->foundOption(OPTION_SCN_DIR)) {

            scn_dir = (char *) cmds->optionValue(OPTION_SCN_DIR).c_str();

        } else {

            fprintf(stderr, "no scenario directory path specified. use "\
                "option -h for help.\n");

            return -1;
        }*/
    }

    if (result == ArgvParser::ParserHelpRequested) {
        return -1;
    }

    // ************************************************************************
    // 2) parse the .scn (scenario) file to extract and/or compute the 
    //    following parameters:
    //      -# LEVEL_DEPTH : nr. of levels
    //      -# LEVEL_BREADTH : nr. of networks/domains per level
    //      -# CONTENT_SOURCE_DIST : nr. of content sources per distance from 
    //          request source (as a function of level)
    //      -# REQUEST_SIZE : nr. of URL elements in the name used to 
    //      -# probability matrix w/ FP probabilities per level
    //      -# o-optimistic values per level
    //      -# latencies and penalties per level
    // ************************************************************************

    DataParser * scn_parser = new DataParser(scn_file);

    // 2.1) nr. of levels
    int level_depth = 0;
    scn_parser->get_int_property_value(LEVEL_DEPTH, level_depth);

    // 2.2) how many networks/domains per level?
    int level_breadth[MAX_LEVEL_DEPTH];
    scn_parser->get_int_property_array(LEVEL_BREADTH, level_breadth);

    // 2.3) nr. of content sources per distance. distances are a 
    // function of the level, e.g. a source which becomes 'visible' at a 
    // level i router will always be (2 * i) + 1 hops away from the origin of
    // the request
    int content_source_dist[MAX_LEVEL_DEPTH];
    scn_parser->get_int_property_array(CONTENT_SOURCE_DIST, content_source_dist);

    // 2.4) request size
    int request_size = 0;
    scn_parser->get_int_property_value(REQUEST_SIZE, request_size);

    // 2.5) distributions of |F\R| per level, a matrix of dimension 
    // LEVEL_DEPTH x REQUEST_SIZE |R|, with the format shown below:
    //
    //           |  1  |  2  | ... | |R| |
    //  ----------------------------------
    //  level 1  | 0.2 | 0.3 | ... | 0.8 | SUM() = 1
    //    ...      ...
    //  level n  | 0.8 | 0.3 | ... | 0.2 | SUM() = 1
    //
    // the values for |F\R| distributions are pre-computed, similarly to 
    // those shown in the qualifier and Hotnets 2015 submissions. the results 
    // are saved in .dat files, in a directory passed as argument.

    // 2.5.1) extract the name of the |F\R| distributions for all levels
    std::string complement_dist[MAX_LEVEL_DEPTH];
    scn_parser->get_string_property_array(COMPLEMENT_DIST, complement_dist);

    // 2.5.2) using a new DataParser, parse the <request_size>.dat file and 
    // get the requested |F\R| distribution for each level
    char dat_file[MAX_FILENAME_SIZE];
    snprintf(dat_file, MAX_FILENAME_SIZE, "%s/%d.dat", dat_dir, request_size);
    DataParser * dat_parser = new DataParser(dat_file);

    double complement_matrix[MAX_LEVEL_DEPTH][MAX_REQUEST_SIZE]; 

    for (int l = 0; l < level_depth; l++) {

        dat_parser->get_double_property_array(complement_dist[l].c_str(), complement_matrix[l]);
    }

    // 2.6) compute the FP probability matrix, of dimension LEVEL_DEPTH X 1 in 
    // the format below:
    //
    //           | P(FP) |
    //  ------------------ 
    //  level 1  |       |
    //  level 2  |       |
    //    ...       ...   
    //  level n  |       |
    //
    // P(FP) is the probability of having the forwarding engine pick a FP entry 
    // at level i. note that this is NOT the same as: 
    //  -# the Bloom Filter FP rate, (1 - exp((m/n) * k))^k, which tells 
    //      us the likelihood of a particular entry to trigger a FP 
    //      match
    // -# the table-wise FP rate, which tells us how likely it is to have 
    //      AT LEAST 1 entry triggering a FP match
    //
    // P(FP) should be calculated according to the Law of Total Probability as:
    //
    // P(FP) = SUM( P(FP | |F\R| = x) * P(|F\R| = x) ), for all x, such that 
    //  1 <= x <= |R|
    //
    // where,
    //  -# P(FP | |F\R| = x) is the Bloom Filter FP rate w/ an |F\R| = x
    //      exponent 
    //  -# P(|F\R| = x) is a point in the distribution of complement sizes, 
    //      passed as a parameter in the .scn file
    //
    double fp_prob[MAX_LEVEL_DEPTH];
    double k = ceil(log(2) * (BF_SIZE / request_size));

    for (int l = 0; l < level_depth; l++) {

        fp_prob[l] = 0.0;

        for (int c = 0; c < request_size; c++) {

            fp_prob[l] += fp_rate((double) BF_SIZE, (double) request_size, k, (double) (c + 1)) * (complement_matrix[l][c] / 100.0);
        }
    }

    // 2.7) o-optimistic transition matrix
    double o_optimistic[MAX_LEVEL_DEPTH][MAX_LEVEL_DEPTH];

    // 2.7.1) initialize everything to 0.0
    for (int i = 0; i < MAX_LEVEL_DEPTH; i++) {
        memset(&o_optimistic[i][0], 0.0, MAX_LEVEL_DEPTH * sizeof(double));
    }

    // 2.7.2) optimistic[0][0] is special, in the sense that we don't specify 
    // the nr. of routers within level 0.
    o_optimistic[0][0] = DEFAULT_O_OPTIMISTIC;

    // 2.7.3) we now build the rest of the o-optimistic matrix, which will give 
    // us the o-optimistic value for a transition between levels. the possible 
    // tranistions are: 
    //  1) i -> i
    //  2) i -> i - 1
    //
    // note that one only goes up when a N decision occurs, and so the 
    // o-optimistic value doesn't come into play for i -> i + 1 transitions

    // case 1 : fill o_optimistic[i][i]. take into account the content sources 
    // which are accessible via level i (minus those which were accessible via 
    // the previous level i - 1 : we assume we never go back), and divide this 
    // by the breadth of the level i - 1 (this approximates the number of 
    // interfaces at level i)
    for (int r = 1; r < level_depth; r++) {
        o_optimistic[r][r] = abs(min((double) (content_source_dist[r] - content_source_dist[r - 1]) / ((double) level_breadth[r - 1] - 1.0), 1.0));
    }

    // case 2 : fill o_optimistic[i][i - 1]. since we're going down from level 
    // i, this is simply 1 / level_breadth[r - 1]
    for (int r = (level_depth - 1); r - 1 >= 0; r--) {

        o_optimistic[r][r - 1] = 1.0 / ((double) level_breadth[r - 1]);
    }

    // 2.7.4) print the transition matrix (later)

    // 2.8) latency values per level, in multiples of some time unit.
    double latencies[MAX_LEVEL_DEPTH];
    scn_parser->get_double_property_array(LATENCIES, latencies);

    // 2.9) the penalty array. our goal is to calculate the relaying 
    // penalty - PENALTY(m) - for a request, given the following conditions: 
    //  -# the request has crossed a level as high as level m
    //  -# we know the distribution of sources per level for the scenario 
    //     being evaluated
    //  -# we know the avg. latencies of forwarding between routers at level l 
    //
    // this can be pre-calculated 'a priori'. how? 
    //
    //  1) for each max. level m, we calculate PENALTY(m, l), i.e. 
    //     the penalties of relaying a request to hypothetical sources on each 
    //     of the other levels l. 
    //  2) then we take the min(PENALTY(m, l)), for all l, and use it as 
    //     PENALTY(m).
    //
    // FIXME: i think i'm making this too complicated
    double penalties[MAX_LEVEL_DEPTH];

    // Q1) say that we know that max. level is m. how do we calculate 
    //     PENALTY(m, l)? we have 3 (m, l) combinations:
    //
    //  -# if l = m : we can have 2 situations: 
    //
    //     1) there may be a source right next to us, reachable with a cost of 
    //        3*L[0]. what is the probability of that happening? it is 
    //     
    //        P(min_cost) = min(content_source_dist[m] / level_breadth[m], 1.0)
    //     
    //     2) other sources may exist which are reachable through level m. in 
    //        this case the cost is L[m] + SUM(2*L[m - i]) + 2*L[0], 
    //        with 1 <= i <= m. the prob. of this happening is 1 - P(min_cost)
    //        
    //        the minimal penalty will *ALWAYS* come from this case, if 
    //        applicable.
    //
    //  -# if l < m : you want to go back to level 0, and through level m.
    //     this always implies a cost of L[m] + SUM(2*L[m - i]) + 2*L[0], 
    //     with 1 <= i <= m. this case provides intermediate penalties.
    //
    //  -# if l > m : you actually have to go through an upper level to 
    //     relay the request. in that case, the cost is L[l] + SUM(2 * L[l - i]) + 2*L[0], 
    //     with 1 <= i <= l. this is the last resort, as it always 
    //     produces the larger penalties.
    // 

    // Q2) how do we know if each one of the 3 cases is applicable 
    // or not?
    // 
    // easy: if there are sources 'visible' at some distance/level l, 
    // we calculate the PENALTY(m, l). 

    // here's a bunch of aux variables, used to keep candidate penalty values 
    // for PENALTY(m), probabilities, etc. 
    double penalty_candidate = DBL_MAX;
    double penalty_min = DBL_MAX;
    double penalty_base = DBL_MAX;

    double prob_penalty_min = 0.0;
    int possible_source_slots = 1;

    for (int m = 0; m < level_depth; m++) {

        // initialize the candidate penalties
        penalty_candidate = DBL_MAX;

        // 
        penalty_base = 2.0 * latencies[0] + latencies[m];
        possible_source_slots = 1;

        for (int l = 0; l < m; l++) {

            penalty_base += 2.0 * latencies[l];
            possible_source_slots = possible_source_slots * level_breadth[l];
        }

        for (int l = 0; l < level_depth; l++) {

            if (content_source_dist[l] < 1)
                continue;

            if (l == m) {

                // prob of finding a source within min_penalty
                prob_penalty_min = min((double) ((double) content_source_dist[m] / (double) possible_source_slots), 1.0);

                // the penalty candidate for this round
                penalty_candidate = (prob_penalty_min * (3.0 * latencies[0])) + ((1.0 - prob_penalty_min) * penalty_base);

            } else if (l > m) {

                penalty_candidate = 2.0 * latencies[0] + latencies[l];

                for (int a = 0; a < l; a++)
                    penalty_candidate += 2.0 * latencies[a];

            } else {

                penalty_candidate = penalty_base;
            }

            // if a new min has been found, update it
            if (penalty_candidate < penalty_min)
                penalty_min = penalty_candidate;
        }

        penalties[m] = penalty_min;
    }

    // 2.9) print all parameters

    printf("\n\n*** PARAMETERS ***\n");

    // 2.9.1) initial row & midrule
    printf("\n[LEVEL]   : |");

    for (int i = 0; i < level_depth; i++)
        printf(" %-8d|", i + 1);

    printf("\n-------------");

    for (int i = 0; i < level_depth; i++)
        printf("----------");

    // 2.9.2) penalty row per level
    printf("\n[SOURCES] : |");

    for (int i = 0; i < level_depth; i++)
        printf(" %-.8d|", content_source_dist[i]);

    // 2.9.3) penalty row per level
    printf("\n[BREADTH] : |");

    for (int i = 0; i < level_depth; i++)
        printf(" %-.8d|", level_breadth[i]);

    // 2.9.4) penalty row per level
    printf("\n[PENALTY] : |");

    for (int i = 0; i < level_depth; i++)
        printf(" %-.6f|", penalties[i]);

    // 2.9.6) o-opt row per level
    printf("\n[FP-PROB] : |");

    for (int i = 0; i < level_depth; i++)
        printf(" %-.2E|", fp_prob[i]);

    // 2.9.5) o-optimistic matrix
    printf("\n\n[O-OPTIMISTIC MATRIX] : \n");
    printf("\n   |");

    for (int i = 0; i < level_depth; i++)
        printf(" %-9d|", i + 1);

    printf("\n-------");

    for (int i = 0; i < level_depth; i++)
        printf("----------");    

    for (int r = 0; r < level_depth; r++) {

        printf("\n %d |", r + 1);

        for (int c = 0; c < level_depth; c++) {
            printf(" %-.3E|", o_optimistic[r][c]);
        }
    }

    printf("\n\n");

    // ************************************************************************
    // 3) build the n-ary tree representing all possible forwarding decisions 
    //    and associated outcomes in the path of a request towards a content 
    //    source, for the scenario(s) described in the .scn file(s)
    // ************************************************************************
    tree<Node *> decision_tree;
    tree<Node *>::iterator root;

    // initialize the probability tree using begin() and add a root 
    // (dummy) node with insert(). after this, we use append_child() all the 
    // way.
    Node * root_node = new Node(0, 0, 0, 1.0, latencies[0], Node::O_NODE);
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

    double next_node_prob = 0.0;
    double path_latency = 0.0;
    double penalty_latency = 0.0;

    // start the .dot file for rendering the probability tree
    Graph * graphviz_graph = new Graph(graphviz_file);

    // depth and breadth are sort of (x, y) coordinates for nodes in a tree. 
    // these will be used to generate node names for the .dot file. the origin 
    // of the tree has coordinates (0, 0).
    int depth = 0;
    int breadth = 0;

    // note the limit condition : on the longest possible path (i.e. going 
    // through the max. level) we make 2 * level_depth decisions
    for (depth = 0; depth < ((2 * level_depth)); depth++) {

        // // just a dummy guard
        // if (depth > 3)
        //     break;

        // update the iterators
        depth_itr = decision_tree.begin_fixed(root, depth);

        breadth = 0;

        // cycle through the previous decisions, append children
        while (decision_tree.is_valid(depth_itr)) {

            curr_node_level = (*depth_itr)->get_next_level();
            curr_node_max_level = (*depth_itr)->get_max_level();

            // if this is a 'dead end', don't append children to it
            if (curr_node_level == END_OF_PATH) {

                printf("skipping children of node : %s(%d, %d)\n", (*depth_itr)->to_string(), depth, breadth);
                
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
                if (content_source_dist[curr_node_level] > 0) {

                    // if there are 'visible' sources, then there's no way we 
                    // can get an N decision. mark it as END_OF_PATH.
                    n_node->set_prob_val(prev_node_prob * 0.0);
                    n_node->set_next_level(END_OF_PATH);

                    // we may then follow a C or I decision. to determine the 
                    // probability of each case, we 'split' them according to 
                    // the o_optimistic value for curr_node_level. 
                    i_node->set_prob_val(prev_node_prob * (1.0 - o_optimistic[curr_node_level][i_node->get_next_level()]) * 1.0);
                    c_node->set_prob_val(prev_node_prob * o_optimistic[curr_node_level][c_node->get_next_level()] * 1.0);

                    // add nodes to the .dot file for graph rendering
                    graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, (1.0 - o_optimistic[curr_node_level][i_node->get_next_level()]) * 1.0);
                    graphviz_graph->add_node((*depth_itr), depth, c_node, depth + 1, breadth, o_optimistic[curr_node_level][c_node->get_next_level()] * 1.0);
                    //graphviz_graph->add_node((*depth_itr), depth, n_node, depth + 1, breadth, 0.0);

                    Node * nodes[] = {i_node, c_node};
                    graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));

                } else {

                    // if there *NO* 'visible' sources, then there's no chance 
                    // of making a C decision: all FPs will be (I)ncorrect. 
                    // mark it as END_OF_PATH.
                    c_node->set_prob_val(prev_node_prob * 0.0);
                    c_node->set_next_level(END_OF_PATH);

                    // we may then follow an I or N decision, which basically 
                    // comes down to P(FP) or (1 - P(FP)).
                    n_node->set_prob_val(prev_node_prob * (1.0 - fp_prob[curr_node_level]));
                    i_node->set_prob_val(prev_node_prob * fp_prob[curr_node_level]);

                    // add nodes to the .dot file for graph rendering
                    graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, fp_prob[curr_node_level]);
                    //graphviz_graph->add_node((*depth_itr), depth, c_node, depth + 1, breadth, 0.0);
                    graphviz_graph->add_node((*depth_itr), depth, n_node, depth + 1, breadth, (1.0 - fp_prob[curr_node_level]));

                    Node * nodes[] = {i_node, n_node};
                    graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
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

                // Q1) if the previous decision was I, we will now go down 
                // a level
                int next_node_level = curr_node_level - 1;

                printf(" *** %s's level is now %d\n", (*depth_itr)->get_uid().c_str(), next_node_level);

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_level == END_OF_PATH) {

                    // penalty due to relaying
                    path_latency = latencies[0];
                    penalty_latency = penalties[curr_node_max_level];

                    // we did not change levels, and so fwd tables will remain 
                    // the same. we'll make the same decision again.
                    next_node_prob = 1.0;

                } else {

                    // the next component of path latency will be that of 
                    // the level below: notice that we're essentially forwarding 
                    // between hops within the same level (in this case, 
                    // next_node_level)
                    path_latency = latencies[next_node_level];
                    penalty_latency = 0.0;

                    next_node_prob = fp_prob[next_node_level];
                }   

                i_node->set_next_level(next_node_level);

                // Q2) the P(FP) we use is that of the next level
                n_node->set_prob_val(prev_node_prob * (1.0 - next_node_prob));
                i_node->set_prob_val(prev_node_prob * next_node_prob);

                // Q3) regarding latencies:
                //  I node : use the values calculated above
                i_node->set_latency_val(prev_node_latency + path_latency + penalty_latency);
                //  N node : the policy is to relay when a I > N transition 
                //           happens, so simply add this level's penalty
                n_node->set_latency_val(prev_node_latency + penalties[curr_node_max_level]);

                // add nodes to the .dot file for graph rendering
                graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, next_node_prob);
                //graphviz_graph->add_node((*depth_itr), depth, c_node, depth + 1, breadth, 0.0);
                graphviz_graph->add_node((*depth_itr), depth, n_node, depth + 1, breadth, (1.0 - next_node_prob));

                Node * nodes[] = {i_node, n_node};
                graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));

            } else {

                // if the previous decision was C, we can either make another 
                // C or an I decision. mark the N node as END_OF_PATH.
                n_node->set_prob_val(prev_node_prob * 0.0);
                n_node->set_next_level(END_OF_PATH);

                // Q1) if the previous decision was C, we will also go down 
                // a level
                int next_node_level = curr_node_level - 1;

                // check if the request is about to be delivered to a content 
                // source
                if (next_node_level == END_OF_PATH) {

                    // penalty due to relaying
                    penalty_latency = penalties[curr_node_max_level];

                    next_node_prob = o_optimistic[0][0];

                } else {

                    penalty_latency = 0.0;
                    next_node_prob = o_optimistic[curr_node_level][next_node_level];
                }

                i_node->set_next_level(next_node_level);
                c_node->set_next_level(next_node_level);

                // Q2) probabilities 
                i_node->set_prob_val(prev_node_prob * (1.0 - next_node_prob) * 1.0);
                c_node->set_prob_val(prev_node_prob * next_node_prob * 1.0);

                // Q3) latencies
                i_node->set_latency_val(prev_node_latency + latencies[next_node_level] + penalty_latency);
                c_node->set_latency_val(prev_node_latency + latencies[next_node_level]);

                // add nodes to the .dot file for graph rendering
                graphviz_graph->add_node((*depth_itr), depth, i_node, depth + 1, breadth, (1.0 - next_node_prob) * 1.0);
                graphviz_graph->add_node((*depth_itr), depth, c_node, depth + 1, breadth, next_node_prob * 1.0);
                //graphviz_graph->add_node((*depth_itr), depth, n_node, depth + 1, breadth, 0.0);

                Node * nodes[] = {i_node, c_node};
                graphviz_graph->align_nodes(nodes, sizeof(nodes) / sizeof(Node *));
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
    graphviz_graph->terminate();

    // the decision tree is finished. let's iterate through the leafs to learn 
    // a bit about latencies for this scenario.
    tree<Node *>::leaf_iterator leaf_itr = decision_tree.begin_leaf();
    tree<Node *>::leaf_iterator end_leaf_itr = decision_tree.end_leaf();

    printf("\n*** FINAL RESULTS ***\n\n");

    double checksum = 0.0;

    while (leaf_itr != end_leaf_itr) {

        if ((*leaf_itr)->get_prob_val() > 0.0) {

            printf("[DEPTH %d] : %s\n", 
                decision_tree.depth(leaf_itr), 
                (*leaf_itr)->to_string());

            checksum += (*leaf_itr)->get_prob_val();
        }

        ++leaf_itr;
    }

    printf("\n[CHECKSUM : %-.8E]\n", checksum);

    return 0;
}
