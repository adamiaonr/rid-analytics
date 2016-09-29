#include <math.h>
#include <cfloat>
#include <algorithm>
#include <string>
#include <fstream>  // reading .nw files
#include <map>
#include <list>

#include "graph.h"
#include "rid-analytics.h"
#include "dataparser.h"
#include "pugixml.hpp"

RID_Analytics::RID_Analytics(
    std::string nw_filename,
    uint8_t request_size,
    uint16_t bf_size) {

    // initialize the scenario parameters set at runtime (e.g. BF size, request 
    // size, TP sizes, F\R distributions)
    this->request_size = request_size;
    this->f_max = request_size;
    this->bf_size = bf_size;
    this->f_r_distribution = NULL;

    // ... and build the network
    this->build_network(nw_filename);
}

RID_Analytics::~RID_Analytics() {

    // erase the path state tree
    tree<Path_State *>::post_order_iterator post_itr = this->path_state_tree.begin_post();
    tree<Path_State *>::post_order_iterator end_post_itr = this->path_state_tree.end_post();

    while (post_itr != end_post_itr) {

        // delete the Path_State object associated with the tree node (not the 
        // node itself)
        delete (*post_itr);

        ++post_itr;
    }

    // erase the routers in the topology. use the router map for that (see the 
    // advantage of using modern C++ features you dummy?)
    for(RID_RouterMap::iterator rit = this->routers.begin(); 
        rit != this->routers.end(); 
        ++rit) {

        delete (*rit).second;
    }
}

int RID_Analytics::get_tp_distance_rec(RID_Router * from_router) {

    int tp_distance = 0;

    // we first check if a TP is present at the cache of from_router
    if (this->tp_sizes[from_router->get_id()][IFACE_LOCAL] > 0)
        return 1;

    // check if it is in a router below from_router
    RID_Router * adj_router = NULL;

    for (int iface = (IFACE_LOCAL + 1); iface < from_router->get_iface_num(); iface++) {

        adj_router = from_router->get_next_hop(iface).router;

        // if the height of the next router is less than or equal to 
        // the current, abort.
        if (adj_router->get_height() < from_router->get_height())
            continue;

        if (adj_router != NULL) {

            if ((tp_distance = get_tp_distance_rec(adj_router)) > 0) {

                return tp_distance + 1;
            }
        }
    }

    return 0;
}

int RID_Analytics::get_tp_distance(RID_Router * from_router) {

    // find the offset of this router's TP size info in the true positive 
    // array
    int height = from_router->get_height();

    // we first check if a TP is present at the cache of from_router
    if (this->tp_sizes[from_router->get_id()][IFACE_LOCAL] > 0)
        return 1;

    // keep a list of visited routers
    std::vector<std::string> visited;

    // if not, check if it is below this router or below its parents. 
    // we start by checking if the true positive lies directly below this 
    // router. if not, we go up one level, and check if it is on another 
    // subtree.
    int tp_distance = 0, go_up = 0;
    RID_Router * parent_router = from_router;
    RID_Router * adj_router = NULL;

    for (go_up = 0; (height - go_up) >= 0; go_up++) {

        // if no more parent routers, abort and return -1 (no TPs in the 
        // topology)
        if (parent_router->get_height() == 0)
            break;

        for (int iface = (IFACE_LOCAL + 1); iface < parent_router->get_iface_num(); iface++) {

            adj_router = parent_router->get_next_hop(iface).router;

            // if the height of the next router is less than the current, abort. 
            // notice that by setting iface = IFACE_LOCAL + 1 we avoid looking 
            // into the local cache and by '<' instead of '<=' we consider 
            // peering relationships.
            if (adj_router->get_height() < parent_router->get_height())
                continue;

            // we've already been on this child, don't go down that road again
            if (std::find(visited.begin(), visited.end(), adj_router->get_id()) != visited.end())
                continue;

            if (adj_router != NULL) {
                if ((tp_distance = get_tp_distance_rec(adj_router)) > 0) {
                    return go_up + tp_distance + 1;
                }
            }
        }

        // looking at this router's subtree didn't yield results. add it to 
        // the list of visited subtrees.
        visited.push_back(parent_router->get_id());
        visited.push_back(adj_router->get_id());
        // update parent_router to its parent
        for (int iface = (IFACE_LOCAL + 1); iface < parent_router->get_iface_num(); iface++) {
            if ((adj_router = parent_router->get_next_hop(iface).router)->get_height() < parent_router->get_height()) {
                parent_router = adj_router;
                break;
            }
        }

        tp_distance = 0;
    }

    // this should never happen (?) this means there are no true positives in 
    // the topology
    return -1;
}

int RID_Analytics::build_network(std::string nw_filename) {

    // read .scn file, structured as an xml file
    pugi::xml_document nw_doc;
    if (!(nw_doc.load_file(nw_filename.c_str()))) {

        std::cerr << "RID_Analytics::get_network() : [ERROR] .scn file (" 
            << nw_filename << ") not found." << std::endl;

        return -1;
    }

    pugi::xml_node trees = nw_doc.child("trees");

    // ***
    // *** cycle through each access tree
    // ***
    for (pugi::xml_node tree = trees.child("tree"); 
        tree; 
        tree = tree.next_sibling("tree")) {

        int tree_index = tree.attribute("index").as_uint();
        this->access_tree_height = tree.attribute("height").as_uint();

        // ***
        // *** extract table sizes, per tier of the access tree
        // ***
        std::cout << "RID_Analytics::build_network() : [INFO] table sizes: ";
        // table sizes saved in a std::map
        std::map<int, __float080> table_sizes;
        for (pugi::xml_node table_size = tree.child("fwd_table_sizes").child("tier"); 
            table_size; 
            table_size = table_size.next_sibling("tier")) {  

            table_sizes[table_size.attribute("tier").as_uint()] = table_size.text().as_uint();

            std::cout << "\n\tTIER = " << table_size.attribute("tier").as_uint()
                << ", SIZE = " << table_size.text().as_uint() 
                    << " (" << table_sizes[table_size.attribute("tier").as_uint()] << ")";
        }
        std::cout << std::endl;

        // ***
        // *** start building the topology by reading the info for each router
        // ***
        for (pugi::xml_node router = tree.child("topology").child("router"); 
            router; 
            router = router.next_sibling("router")) {

            int tier = router.attribute("tier").as_uint();
            int tier_index = router.attribute("index").as_uint();
            // the router id follows the format <tree index>.<tier>.<tier index>
            std::string router_id = std::string(tree.attribute("index").value())
                + std::string(".")
                + std::string(router.attribute("id").value());

            // ***
            // *** extract iface nr. & entry distr. per iface
            // ***
            std::cout << "RID_Analytics::build_network() : [INFO] \% of entries per iface: ";
            // extract the nr. of ifaces by counting the number of fwd_dist xml nodes
            int iface_num = 0;
            // extract iface <-> map
            std::map<int, __float080> fwd_dist;
            for (pugi::xml_node dist = router.child("fwd_dist"); 
                dist; 
                dist = dist.next_sibling("fwd_dist")) {

                fwd_dist[(int) dist.attribute("iface").as_int()] = ((__float080) dist.text().as_double());
                // increment the number of ifaces here
                iface_num++;

                std::cout << "\n\tIFACE = " << dist.attribute("iface").as_int()
                    << ", DISTR. = " << (__float080) dist.text().as_double() 
                        << " (" << fwd_dist[(int) dist.attribute("iface").as_int()] << ")";
            }
            std::cout << std::endl;

            // ***
            // *** create or initialize new RID router
            // ***
            RID_Router * new_router = NULL;
            RID_RouterMap::iterator rit;
            if ((rit = this->routers.find(router_id)) != this->routers.end()) {

                new_router = (*rit).second;
                new_router->init(
                    tier, tier_index,           // (tier, left-to-right index) coordinates
                    table_sizes[tier],          // nr. of table entries for tier
                    iface_num,                  // nr. of ifaces equals nr. of adjacencies
                    this->f_max,                // max. request & forwarding entry size
                    this->bf_size);             // BF size (of both requests and entries, in bit)

            } else {

                new_router = 
                    new RID_Router(
                        tree_index,
                        this->access_tree_height,
                        tier, tier_index,
                        table_sizes[tier],
                        iface_num,
                        this->f_max,
                        this->bf_size);

                // add the router to the router map
                this->routers[router_id] = new_router;
            }

            // ***
            // *** initialize router ifaces & add fwd entries
            // ***
            for (pugi::xml_node link = router.child("link"); 
                link; 
                link = link.next_sibling("link")) {

                // ***
                // *** extract fwd entry size distributions
                // ***
                std::map<int, __float080> size_dist;
                for (pugi::xml_node s_dist = link.child("fwd_size_dist"); 
                    s_dist; 
                    s_dist = s_dist.next_sibling("fwd_size_dist")) {

                    size_dist[s_dist.attribute("size").as_uint()] = s_dist.text().as_double();
                }

                RID_Router * adjacent_router = NULL;
                RID_RouterMap::iterator arit;
                // if adjacent router is already in the router map, fetch its 
                // pointer. otherwise, create a new RID_Router object.
                std::string adjacent_router_id = std::string(tree.attribute("index").value())
                    + std::string(".") 
                    + std::string(link.attribute("rrouter").value());

                if ((arit = this->routers.find(adjacent_router_id)) != this->routers.end()) {

                    adjacent_router = (*arit).second;

                } else {

                    adjacent_router = new RID_Router(tree_index, this->access_tree_height);
                    // add adj. router to the router map
                    this->routers[adjacent_router_id] = adjacent_router;
                }

                new_router->add_fwd_table_entry(
                    link.attribute("local").as_uint(),
                    fwd_dist[link.attribute("local").as_uint()],
                    size_dist,
                    adjacent_router, link.attribute("remote").as_uint());
            }
        }
    }

    return 0;
}

int RID_Analytics::run_rec(
    RID_Router * router,
    uint8_t ingress_iface,
    tree<Path_State *>::iterator prev_path_state_itr) {

    // local vars for quick access to router attributes
    uint8_t iface_num = router->get_iface_num();
    uint8_t height = router->get_height();
    uint8_t width = router->get_width();

    __float080 event_prob = 0.0;
    __float080 ingress_prob = (*prev_path_state_itr)->get_ingress_iface_prob();
    __float080 * ingress_ptree_prob = (*prev_path_state_itr)->get_ingress_ptree_prob();

    // next hop information
    RID_Router::nw_address next_hop;
    tree<Path_State *>::iterator path_state_itr;

    // run the forward() operation on this router : this will determine the 
    // following info:
    //  1) interface event probabilities
    //  2) egress size probabilities
    //
    // we then use this info to determine the path state
    printf("\nRID_Analytics::run_rec() : FORWARD() BY ROUTER[%d][%d]\n", height, width);

    router->forward(
        this->request_size, 
        ingress_iface,
        this->tp_sizes[router->get_id()],
        ingress_prob,
        ingress_ptree_prob,
        this->f_r_distribution);

    // we now cycle through all possible interface events and set the possible 
    // path states and probabilities
    for (int event = EVENT_NO_MATCHES; event < EVENT_NUM; event++) {

        // we differentiate the establishment of states by iface events:
        //  1) if EVENT_NO_MATCHES, EVENT_MULTIPLE_IFACE_MATCHES 
        //     or EVENT_LOCAL_IFACE_MATCH the path ends
        //  2) if EVENT_SINGLE_IFACE_MATCH, the path continues, and so we have 
        //     additional work : setting ingress 'prefix tree' probabilities, 
        //     determining the intermediate state (TP or FP), calling run_rec() 
        //     again, etc.
        Path_State * path_state = NULL;

        if (event < EVENT_SINGLE_IFACE_MATCH) {

            path_state = new Path_State(router, this->request_size);
            path_state_itr = this->path_state_tree.append_child(prev_path_state_itr, path_state);

            // all events < EVENT_SINGLE_IFACE_MATCH are 'ends of the line' 
            // in a way.
            path_state->set_final_prob(router->get_iface_events_prob(event));
            path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

            if (event > EVENT_NO_MATCHES) {

                // as long as the event is EVENT_MULTIPLE_IFACE_MATCHES or
                // EVENT_LOCAL_IFACE_MATCH, the path ends here...
                // FIXME: this is bound to change.
                path_state->set_eop();

                if (event == EVENT_MULTIPLE_IFACE_MATCHES) {

                    path_state->set_outcome(OUTCOME_FALLBACK_DELIVERY);
                    path_state->set_path_length(
                        (*prev_path_state_itr)->get_path_length() + 1);
                        // + get_tp_distance(path_state->get_router()));

                    printf("\nRID_Analytics::run_rec() : Event 'MIS' at ROUTER[%d][%d] : latency = %d\n", 
                        height, width, path_state->get_path_length());

                } else if (event == EVENT_LOCAL_IFACE_MATCH) {

                    // if we hit a cache, and there's a TP at iface local, 
                    // we get a correct delivery, otherwise an incorrect delivery
                    if (this->tp_sizes[router->get_id()][IFACE_LOCAL] > 0) {

                        path_state->set_outcome(OUTCOME_CORRECT_DELIVERY);

                    } else {

                        path_state->set_outcome(OUTCOME_INCORRECT_DELIVERY);
                        // this is the fancy TP distance thing (in fact, not 
                        // that fancy and/or useful)
                        path_state->set_path_length(
                            (*prev_path_state_itr)->get_path_length() + 1);
                            // + get_tp_distance(path_state->get_router()));
                    }
                }

            } else {

                // EVENT_NO_MATCHES leads to the relay of the packet
                if (ingress_iface != IFACE_UPSTREAM 
                    && (router->get_height() > 0 && router->get_width() > 0)) {
                    path_state->set_final_prob(0.0);
                }

                path_state->set_outcome(OUTCOME_FALLBACK_RELAY);
                path_state->set_path_length(
                    (*prev_path_state_itr)->get_path_length() + 1);
                    // + get_tp_distance(path_state->get_router()));
            }

        } else {

            // if EVENT_SINGLE_IFACE_MATCH, the request paths continue to be 
            // built. therefore, we cycle through all egress interfaces
            for (int iface = IFACE_LOCAL + 1; iface < iface_num; iface++) {

                // create a new PathState object for each egress iface
                path_state = new Path_State(router, this->request_size);
                path_state_itr = this->path_state_tree.append_child(prev_path_state_itr, path_state);

                // increment the path length by +1 hop
                path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

                // now, we deal with probabilities:
                //  -# ingress_ptree_probs : the prob of having the packet 
                //     bound to a prefix tree of size s (no TP info)
                //
                //  -# ingress_iface_probs : the prob of having a packet 
                //     flow through iface i (takes TPs into account)
                //
                path_state->set_ingress_ptree_prob(router->get_egress_ptree_prob(iface), this->f_max);

                event_prob = router->get_egress_iface_prob(iface);
                path_state->set_ingress_iface_prob(event_prob);
                path_state->set_final_prob(event_prob);

                // for (int f = 0; f < this->f_max + 1; f++) {
                //     event_prob = router->get_egress_iface_prob(iface, f);
                //     cumulative_event_prob += event_prob;
                //     path_state->set_ingress_ptree_prob(f, router->get_egress_ptree_prob(iface, f));
                //     path_state->set_ingress_prob(f, cumulative_event_prob);
                // }

                // determine the correctness of the forwarding decision
                if (this->tp_sizes[router->get_id()][iface] > 0) {

                    // if router has a TP on iface, set the even as 
                    // an intermediate TP, i.e. "we're on the right path"
                    path_state->set_outcome(OUTCOME_INTERMEDIATE_TP); 

                } else {

                    // else, an 'unlucky' FP happened: we're on the wrong path 
                    // now (a FP match occurred, in an iface which does not 
                    // have TP entries)
                    path_state->set_outcome(OUTCOME_INTERMEDIATE_FP);
                }

                // finally, we continue with the request path
                // get the next hop router by iface
                next_hop = router->get_next_hop(iface);

                // this can still happen (or can it?)
                // FIXME: this subtle condition is the one that actually 
                // terminates the run. this could be the source of trouble.
                if (next_hop.router == router || iface == ingress_iface)
                    continue;

                // if the RID packet came from upwards or a peer router, 
                // make sure to ONLY forward downwards.
                if ((router->get_next_hop(ingress_iface).router->get_height() <= router->get_height())
                        && (router->get_id() != "0.3.0")) {

                    if(next_hop.router->get_height() <= router->get_height())
                        continue;
                }

                std::cout << "RID_Analytics::run_rec() : [INFO] on router[" << router->get_id() 
                    << "], forwarding to router[" << next_hop.router->get_id() 
                    << "]" << std::endl;

                run_rec(next_hop.router, next_hop.iface, path_state_itr);
            }
        }
    }

    return 0;
}

int RID_Analytics::run(std::string nw_filename) {

    // extract the TP size map and |F\R| distr. from nw_filename
    // read .scn file, structured as an xml file
    pugi::xml_document nw_doc;
    if (!(nw_doc.load_file(nw_filename.c_str()))) {

        std::cerr << "RID_Analytics::run() : [ERROR] .scn file (" 
            << nw_filename << ") not found." << std::endl;

        return -1;
    }

    pugi::xml_node trees = nw_doc.child("trees");
    for (pugi::xml_node tree = trees.child("tree"); 
        tree; 
        tree = tree.next_sibling("tree")) {

        // tp sizes per router
        for (pugi::xml_node router = tree.child("topology").child("router"); 
            router; 
            router = router.next_sibling("router")) {

            std::string router_id = std::string(tree.attribute("index").value()) 
                + std::string(".") 
                + std::string(router.attribute("id").value());

            int * _tp_sizes = (int *) calloc(routers[router_id]->get_iface_num(), sizeof(int));
            for (pugi::xml_node tp = router.child("tp"); 
                tp; 
                tp = tp.next_sibling("tp")) {

                _tp_sizes[tp.attribute("iface").as_uint()] = tp.text().as_uint();
            }

            // we keep tp sizes in a map router_id <-> tp sizes (per iface)
            this->tp_sizes[router_id] = _tp_sizes;
        }

        this->f_r_distribution = (__float080 *) calloc(this->f_max, sizeof(__float080));
        for (pugi::xml_node f_r_dist = tree.child("topology").child("f_r_dist"); 
            f_r_dist; 
            f_r_dist = f_r_dist.next_sibling("f_r_dist")) {

            this->f_r_distribution[f_r_dist.attribute("diff").as_uint()] = f_r_dist.text().as_uint();
        }
    }

    tree<Path_State *>::iterator path_state_tree_itr;
    path_state_tree_itr = this->path_state_tree.begin();

    // FIXME: not sure if it is ok to set the rid_router arg on Path_State() 
    // as NULL
    Path_State * initial_state = new Path_State(NULL, request_size);
    initial_state->set_path_length(0);

    path_state_tree_itr = this->path_state_tree.insert(
                                                    path_state_tree_itr, 
                                                    initial_state);

    // set the initial ingress probabilities: ptree_size = 0 gets 1.0 
    // probability, i.e. initially the request isn't bound to any prefix tree
    initial_state->set_ingress_ptree_prob(0, 1.0);
    initial_state->set_ingress_iface_prob(1.0);

    // all other prefix tree sizes get 0.0 probability
    for (uint8_t ptree_size = 1; ptree_size <= this->request_size; ptree_size++)
        initial_state->set_ingress_ptree_prob(ptree_size, 0.0);

    // ingress iface is the IFACE_DOWNSTREAM of upstream router
    run_rec(this->routers["0.3.0"], 2, path_state_tree_itr);

    return 0;
}

int RID_Analytics::view_results(
    uint8_t modes,
    char * output_file_path) {

    // iterate through the path_tree leafs to learn 
    // a bit about latencies for this scenario
    tree<Path_State *>::leaf_iterator _leaf_itr = this->path_state_tree.begin_leaf();
    tree<Path_State *>::leaf_iterator _end_leaf_itr = this->path_state_tree.end_leaf();

    __float080 _checksum = 0.0;

    char * _path_state_str = NULL;

    // print some results for quick checking
    if ((modes & MODE_VERBOSE))
        printf("\n*** RESULTS ***\n\n");

    // if requested, save outcomes in the path given as arg
    FILE * _output_file;

    if ((modes & MODE_SAVE_OUTCOMES)) {

        _output_file = fopen(output_file_path, "w");
    }

    while (_leaf_itr != _end_leaf_itr) {

        if ((modes & MODE_VERBOSE)) {

            _path_state_str = (*_leaf_itr)->to_string();

            printf("[PATH_TREE DEPTH %d] : %s\n", 
                this->path_state_tree.depth(_leaf_itr), 
                _path_state_str);

            free(_path_state_str);
        }

        if ((modes & MODE_SAVE_OUTCOMES)) {

            // append a line in the output file, with the following 
            // format : 
            // <rtr.hght>,<rtr.wdth>,<prob>,<path_length>,<outcome>
            fprintf(_output_file, 
                "%d,%d,%-.5LE,%d,%d\n", 
                (*_leaf_itr)->get_router()->get_height(),
                (*_leaf_itr)->get_router()->get_width(),
                (*_leaf_itr)->get_final_prob(),
                (*_leaf_itr)->get_path_length(),
                (*_leaf_itr)->get_outcome());
        }

        _checksum += (*_leaf_itr)->get_final_prob();
        ++_leaf_itr;
    }

    if ((modes & MODE_VERBOSE)) {

        printf("\n[CHECKSUM : %-.5LE]\n", _checksum);
    }

    return 0;
}
