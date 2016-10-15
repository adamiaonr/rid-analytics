#include <math.h>
#include <cfloat>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>  // reading .nw files
#include <map>
#include <list>

#include "graph.h"
#include "rid-analytics.h"
#include "dataparser.h"
#include "pugixml.hpp"

const char * PATH_TO_ORIGIN[] = { 
    "0.3.0:1", 
    "0.2.0:1", 
    "0.1.0:1", 
    "0.0.0:2", 
    "0.1.1:3", 
    "0.2.3:3", 
    "0.3.7:0", };
const int PATH_TO_ORIGIN_SIZE = 7;

RID_Analytics::RID_Analytics(
    std::string nw_filename,
    uint8_t request_size,
    uint16_t bf_size,
    std::string origin_server,
    int mm_mode,
    int eh_mode) {

    // initialize the scenario parameters set at runtime (e.g. BF size, request 
    // size, TP sizes, F\R distributions)
    this->request_size = request_size;
    this->f_max = request_size;
    this->bf_size = bf_size;
    this->mm_mode = mm_mode;
    this->eh_mode = eh_mode;
    this->f_r_distribution = NULL;

    // ... and build the network
    this->build_network(nw_filename);
    this->read_scn(nw_filename);

    // after building the network, save the origin server pointer
    this->origin_server = this->routers[origin_server];
    // initialize the path to origin, useful for fallbacks
    // FIXME: the path to origin is temporarily hardcoded
    for (int i = 0; i < PATH_TO_ORIGIN_SIZE; i++)
        this->origin_server_path.insert(std::string(PATH_TO_ORIGIN[i]));
}

RID_Analytics::~RID_Analytics() {

    // erase the path state tree
    tree<Path_State *>::iterator dfs_itr = this->path_state_tree.begin();
    // tree<Path_State *>::post_order_iterator post_itr = this->path_state_tree.begin_post();
    // tree<Path_State *>::post_order_iterator end_post_itr = this->path_state_tree.end_post();

    while (dfs_itr != this->path_state_tree.end()) {

        if ((*dfs_itr)->get_router()) {
            std::cout << "RID_Analytics::get_network() : [INFO] deleting {"
                << (*dfs_itr)->get_path_length() << ", "
                << (*dfs_itr)->get_router()->get_id() << ", "
                << (*dfs_itr)->get_path_status() << ", "
                << (*dfs_itr)->get_path_prob() << "}" << std::endl;
        }

        // delete the Path_State object associated with the tree node (not the 
        // node itself)
        delete(*dfs_itr);
        ++dfs_itr;
    }

    // erase the routers in the topology. use the router map for that (see the 
    // advantage of using modern C++ features you dummy?)
    for(RID_RouterMap::iterator rit = this->routers.begin(); 
        rit != this->routers.end(); 
        ++rit) {

        delete (*rit).second;
    }
}

int RID_Analytics::get_origin_distance_rec(RID_Router * from_router) {

    int tp_distance = 0;

    // check if we're there yet
    if (from_router->get_id() == this->origin_server->get_id())
        return 0;

    // check if it is in a router below from_router
    RID_Router * adj_router = NULL;

    for (int iface = (IFACE_LOCAL + 1); iface < from_router->get_iface_num(); iface++) {

        adj_router = from_router->get_next_hop(iface).router;

        // if the height of the next router is less than or equal to 
        // the current, abort.
        if (adj_router->get_height() < from_router->get_height())
            continue;

        if (adj_router->get_id() == from_router->get_id())
            continue;

        if (adj_router != NULL) {

            if ((tp_distance = get_origin_distance_rec(adj_router)) >= 0) {

                return tp_distance + 1;
            }
        }
    }

    return -1;
}

int RID_Analytics::get_origin_distance(RID_Router * from_router) {

    // check if we're there yet
    if (from_router->get_id() == this->origin_server->get_id())
        return 0;

    // if not, check if it is below this router or below its parents. 
    // we start by checking if the true positive lies directly below this 
    // router. if not, we go up one level, and check if it is on another 
    // subtree.
    int tp_distance = 0, go_up = 0, height = from_router->get_height();
    RID_Router * parent_router = from_router;
    RID_Router * adj_router = NULL;
    // keep a list of visited routers
    std::vector<std::string> visited;
    visited.push_back(from_router->get_id());

    for (go_up = 0; (height - go_up) >= 0; go_up++) {

        if ((height - go_up) == 0) {

            std::cout << "RID_Analytics::get_origin_distance() : [INFO] dist. from " 
                << from_router->get_id() << " to " << origin_server->get_id() << " is " << (int) (go_up + 3) << " hops" << std::endl;

            return go_up + 3;
        }

        for (int iface = (IFACE_LOCAL + 1); iface < parent_router->get_iface_num(); iface++) {

            adj_router = parent_router->get_next_hop(iface).router;

            // if the height of the next router is less than the current, skip. 
            // notice that by setting iface = IFACE_LOCAL + 1 we avoid looking 
            // into the local cache and by '<' instead of '<=' we consider 
            // peering relationships.
            if (adj_router->get_height() < parent_router->get_height())
                continue;

            // we've already been on this child, don't go down that road again
            if (std::find(visited.begin(), visited.end(), adj_router->get_id()) != visited.end())
                continue;

            if (adj_router->get_id() == parent_router->get_id())
                continue;

            if (!(parent_router->get_num_entries(iface) > 0))
                continue;

            if (adj_router != NULL) {
                if ((tp_distance = get_origin_distance_rec(adj_router)) >= 0) {

                    std::cout << "RID_Analytics::get_origin_distance() : [INFO] dist. from " 
                        << from_router->get_id() << " to " << origin_server->get_id() << " is " << (int) (go_up + tp_distance + 1) << " hops" << std::endl;

                    return go_up + tp_distance;
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
                    this->bf_size,              // BF size (of both requests and entries, in bit)
                    this->mm_mode);

            } else {

                new_router = 
                    new RID_Router(
                        tree_index,
                        this->access_tree_height,
                        tier, tier_index,
                        table_sizes[tier],
                        iface_num,
                        this->f_max,
                        this->bf_size,
                        this->mm_mode);

                // add the router to the router map
                this->routers[router_id] = new_router;
            }

            if (tier == this->access_tree_height && tier_index > 0)
                new_router->set_leaf();

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

int RID_Analytics::on_path_to_origin(RID_Router * router, int iface) {

    if ((iface == -1)) {

        for (std::set<std::string>::iterator it = this->origin_server_path.begin(); 
            it != this->origin_server_path.end();
            ++it) {

            std::size_t found = (*it).find(router->get_id());
            if (found != std::string::npos)
                return 1;
        }

        return 0;
    }

    std::string the_situation = router->get_id() + std::string(":") + std::to_string(iface);
    std::set<std::string>::iterator it = this->origin_server_path.find(the_situation);
    if (it != this->origin_server_path.end())
        return 1;
    else
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

    __float080 path_prob = 0.0, fallback_carry_prob = 0.0;
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
    std::cout << "\nRID_Analytics::run_rec() : FORWARD() BY ROUTER[" 
        << (int) height << "][" << (int) width << "]" << std::endl;

    if (this->tp_sizes.empty())
        std::cout << "RID_Analytics::run_rec() : [INFO] TP sizes is empty!";
    else {

        std::cout << "RID_Analytics::run_rec() : [INFO] TP sizes keys : " << std::endl;
        RID_TPMap::iterator itr;
        for (itr = this->tp_sizes.begin(); itr != this->tp_sizes.end(); itr++)
            std::cout << "\tTP[" << itr->first << "]" << std::endl;
    }

    // for (int i = 0; i < routers[router->get_id()]->get_iface_num(); i++)
    //     std::cout << "\tTP[" << i << "] = " << (int) this->tp_sizes[router->get_id()][i] << std::endl;       

    router->forward(
        this->request_size, 
        ingress_iface,
        this->tp_sizes[router->get_id()],
        ingress_prob,
        ingress_ptree_prob,
        this->f_r_distribution);

    // we now cycle through all possible interface events and set the possible 
    // path states and probabilities
    for (int event = EVENT_NLM; event < EVENT_NUM; event++) {

        // we differentiate the establishment of states by iface events:
        //
        //  1) if EVENT_SINGLE_IFACE_MATCH, the path continues, and so we have 
        //     additional work : setting FP tree probabilities, 
        //     determining the intermediate state (TP or FP), calling run_rec() 
        //     again, etc.
        //
        //  2) otherwise, the path may or may not end
        //
        if (event < EVENT_SLM) {

            Path_State * path_state = new Path_State(router, this->request_size);
            path_state_itr = this->path_state_tree.append_child(prev_path_state_itr, path_state);

            // set the event & event probability
            path_state->set_event(event, router->get_iface_events_prob(event));
            // all events < EVENT_SINGLE_IFACE_MATCH are 'ends of the line' 
            // in a way.
            path_state->set_path_prob(router->get_iface_events_prob(event));
            path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

            std::cout << "RID_Analytics::run_rec() : [INFO] ADDED NODE TO PATH STATE TREE: "
                << "\n\tEVENT : " << event << " EVENT PROB. : " << path_state->get_event_prob()
                << "\n\tPROB : " << router->get_iface_events_prob(event)
                << "\n\tPATH LATENCY : " << (*prev_path_state_itr)->get_path_length() + 1
                << std::endl;

            if (event > EVENT_NLM) {

                path_state->set_eop();

                if (event == EVENT_MLM) {

                    if (this->mm_mode == MMH_FALLBACK) {

                        path_state->set_path_status(OUTCOME_FALLBACK_DELIVERY);

                        // FIXME: this goes against the on-path resolution of 
                        // fallbacks, but it helps in getting results ready. 
                        // this should be revised in the future.
                        if (on_path_to_origin(router, -1)) {

                            // if the router is on the path to origin, we copy the 
                            // probability of multiple link matches to the egress link 
                            // used to get to the origin server...
                            fallback_carry_prob += router->get_iface_events_prob(event);
                            // ... and set the path length of this state to 0.0
                            // FIXME: this shouldn't be necessary
                            path_state->set_path_length(0.0);

                        } else {

                            // otherwise, we set the path length for this MLM 
                            // event
                            path_state->set_path_length(
                                (*prev_path_state_itr)->get_path_length() + 1);

                            // we add a new state representing the fallback 
                            // delivery, appending it to the MLM event state
                            //  -# the router of the resolution state should 
                            //     be the origin server
                            //  -# append it to 'path_state_itr'
                            //  -# event should be EVENT_LLM 
                            //  -# event prob should be the same as P(MLM)
                            //  -# path status should be 'OUTCOME_FALLBACK_DELIVERY'
                            //  -# path prob should be the same as P(MLM)
                            //  -# path length is set to current length + 1 + distance to origin server
                            Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                            this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                            resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                            resolution_state->set_path_status(OUTCOME_FALLBACK_DELIVERY);
                            resolution_state->set_path_prob(router->get_iface_events_prob(event));
                            resolution_state->set_path_length(
                                (*prev_path_state_itr)->get_path_length() + 1 + get_origin_distance(router));
                        }

                    } else {

                        // if other handling method other than MMH_FALLBACK, then 
                        // the prob of EVENT_MLM will be transfered to STATUS_TP 
                        // or STATUS_FP
                        path_state->set_path_prob(0.0);
                    }

                } else if (event == EVENT_LLM) {

                    // if we hit a cache, and there's a TP at iface local, 
                    // we get a correct delivery, otherwise an incorrect delivery
                    if (this->tp_sizes[router->get_id()][IFACE_LOCAL] > 0) {

                        path_state->set_path_status(OUTCOME_CORRECT_DELIVERY);

                    } else {

                        path_state->set_path_status(OUTCOME_INCORRECT_DELIVERY);
                        path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

                        int penalty = 0;
                        if (this->mm_mode == MMH_FALLBACK) {

                            if (on_path_to_origin(router, -1)) {

                                fallback_carry_prob += router->get_iface_events_prob(event);
                                // ... and set the path length of this state to 0.0
                                // FIXME: this shouldn't be necessary
                                path_state->set_path_length(0.0);

                            } else {

                                penalty = get_origin_distance(router);

                                Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                                this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                                resolution_state->set_path_status(OUTCOME_FALLBACK_DELIVERY);
                                resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                                resolution_state->set_path_prob(router->get_iface_events_prob(event));
                                resolution_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1 + penalty);
                            }

                        } else {

                            // otherwise, the request source must receive 
                            // feedback of the incorrect delivery, in addition 
                            // to the time of having the request 
                            // fetched from the origin server.
                            penalty = (*prev_path_state_itr)->get_path_length() 
                                + get_origin_distance(this->routers["0.3.0"]);

                            Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                            this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                            resolution_state->set_path_status(OUTCOME_CORRECT_DELIVERY);
                            resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                            resolution_state->set_path_prob(router->get_iface_events_prob(event));
                            resolution_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1 + penalty);
                        }
                    }
                }

            } else {

                // if EVENT_NLM happens as the packet is going up the topology, 
                // then it is relayed to the upper tier of the topology. as 
                // such, the probability of this event will transfer to 
                // the SLM event. for this reason, we set P^{E}(NLM) = 0.0
                if (ingress_iface != IFACE_UPSTREAM 
                    && (router->get_height() > 0 && router->get_width() > 0)) {

                    path_state->set_path_prob(0.0);
                    path_state->set_path_status(STATUS_TN);

                } else {

                    int penalty = 0;
                    if (this->mm_mode == MMH_FALLBACK) {

                        if (on_path_to_origin(router, -1)) {

                            fallback_carry_prob += router->get_iface_events_prob(event);
                            // ... and set the path length of this state to 0.0
                            // FIXME: this shouldn't be necessary, now to calculate 
                            // avg. latencies you should only look into correct 
                            // deliveries
                            path_state->set_path_length(0.0);

                        } else {

                            penalty = get_origin_distance(router);

                            Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                            this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                            resolution_state->set_path_status(OUTCOME_FALLBACK_DELIVERY);
                            resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                            resolution_state->set_path_prob(router->get_iface_events_prob(event));
                            resolution_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1 + penalty);
                        }

                    } else {

                        penalty = (*prev_path_state_itr)->get_path_length() 
                                    + get_origin_distance(this->routers["0.3.0"]);

                        path_state->set_path_status(OUTCOME_PACKET_DROP);
                        path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

                        Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                        this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                        resolution_state->set_path_status(OUTCOME_CORRECT_DELIVERY);
                        resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                        resolution_state->set_path_prob(router->get_iface_events_prob(event));
                        resolution_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1 + penalty);
                    }
                }
            }

        } else {

            bool event_added = false;

            // if EVENT_SINGLE_IFACE_MATCH, the request paths continue to be 
            // built. therefore, we cycle through all egress interfaces
            for (int iface = IFACE_LOCAL + 1; iface < iface_num; iface++) {

                // don't add a state for a blocked iface
                std::set<int> blocked_ifaces = router->get_blocked_ifaces();
                std::set<int>::iterator it = blocked_ifaces.find((int) iface);
                if (it != blocked_ifaces.end())
                    continue;

                // create a new PathState object for each egress iface
                Path_State * path_state = new Path_State(router, this->request_size);
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

                path_prob = router->get_egress_iface_prob(iface);

                // FIXME: this only works for the opportunistic caching 
                // scenario

                if (on_path_to_origin(router, iface) && (this->mm_mode == MMH_FALLBACK)) {

                    std::cout << "RID_Analytics::run_rec() : [INFO] router " 
                        << router->get_id() << " is on the path to the origin "
                        << ", adding " << fallback_carry_prob << " to P^E(SLM)[" 
                        << iface << "] = " << path_prob << std::endl;

                    path_prob += fallback_carry_prob;
                }

                path_state->set_path_prob(path_prob);
                path_state->set_ingress_iface_prob(path_prob);

                // set the event & event probability
                if (!event_added) {
                    path_state->set_event(event, router->get_iface_events_prob(event) + fallback_carry_prob);
                    event_added = true;
                } else {
                    path_state->set_event(event, 0.0);
                }

                std::cout << "RID_Analytics::run_rec() : [INFO] ADDED NODE TO PATH STATE TREE: "
                    << "\n\tEVENT : " << event << " EVENT PROB. : " << path_state->get_event_prob()
                    << "\n\tPATH PROB. : " << path_prob
                    << "\n\tPATH LATENCY : " << (*prev_path_state_itr)->get_path_length() + 1
                    << std::endl;

                // determine the correctness of the forwarding decision
                if (this->tp_sizes[router->get_id()][iface] > 0) {

                    // if router has a TP on iface, set the even as 
                    // an intermediate TP, i.e. "we're on the right path"
                    path_state->set_path_status(STATUS_TP); 

                } else {

                    // else, an 'unlucky' FP happened: we're on the wrong path 
                    // now (a FP match occurred, in an iface which does not 
                    // have TP entries)
                    path_state->set_path_status(STATUS_FP);
                }

                // finally, we continue with the request path
                // get the next hop router by iface
                next_hop = router->get_next_hop(iface);

                // this can still happen (or can it?)
                // FIXME: this subtle condition is the one that actually 
                // terminates the run. 
                // FIXME: this could be the source of trouble.
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

int RID_Analytics::read_scn(std::string nw_filename) {

    // extract the TP size map and |F\R| distr. from nw_filename
    // read .scn file, structured as an xml file
    pugi::xml_document nw_doc;
    if (!(nw_doc.load_file(nw_filename.c_str()))) {

        std::cerr << "RID_Analytics::read_scn() : [ERROR] .scn file (" 
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

            int i = (f_r_dist.attribute("diff").as_int() - 1);
            this->f_r_distribution[i] = (__float080) f_r_dist.text().as_double();

            std::cout << "RID_Analytics::read_scn() : [INFO] |F\\R|[" 
                << i << "] = " << this->f_r_distribution[i] << std::endl;
        }
    }

    return 0;
}

int RID_Analytics::run(std::string nw_filename) {

    if (this->tp_sizes.empty())
        std::cout << "RID_Analytics::run() : [INFO] TP sizes is empty!";
    else {

        std::cout << "RID_Analytics::run() : [INFO] TP sizes keys : " << std::endl;
        RID_TPMap::iterator itr;
        for (itr = this->tp_sizes.begin(); itr != this->tp_sizes.end(); itr++)
            std::cout << "\tTP[" << itr->first << "]" << std::endl;
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

std::string get_unix_timestamp() {

    std::stringstream strm;
    strm << time(NULL);

    return strm.str();
}

int RID_Analytics::view_results(
    uint8_t modes,
    std::string output_dir,
    std::string output_label) {

    __float080 checksum = 0.0;
    char * path_state_str = NULL;

    // if verbose mode is on, print results in stdout
    if ((modes & MODE_VERBOSE))
        std::cout << "\n*** RESULTS ***\n" << std::endl;

    // we create x .tsv files, each for different
    //  -# events.tsv       : probabilities of events (LIM, SIM, NIM, MIM) per AS
    //  -# outcomes.tsv     : outcomes (correct/incorrect delivery, etc.) per AS
    //  -# latencies.tsv    : final path latencies (nr. of hops)
    ofstream output_file[2];  
    // file names follow the convention <type>.<label>.<unix-timestamp>.tsv
    std::string output_filename[2] = {"events", "path"};
    for (int i = 0; i < 2; i++) {
        output_filename[i] += std::string(".") + output_label + std::string(".") + get_unix_timestamp() + std::string(".tsv");
        output_file[i].open(output_dir + std::string("/") + output_filename[i]);
    }

    // add column lines to .tsv files
    output_file[FILE_EVENTS] << "AS\tTIER-Y\tTIER-X\tEVENT\tPROB\n";
    output_file[FILE_PATHS] << "AS\tTIER-Y\tTIER-X\tSTATUS\tLATENCY\tPROB\n";

    // gather event stats (DFS iterator traverses the whole path state tree)
    tree<Path_State *>::iterator dfs_itr = this->path_state_tree.begin();
    while(dfs_itr != this->path_state_tree.end()) {

        if ((*dfs_itr)->get_router())
            output_file[FILE_EVENTS] << (*dfs_itr)->get_router()->get_id() << "\t"
                << (int) (*dfs_itr)->get_router()->get_height() << "\t" << (int) (*dfs_itr)->get_router()->get_width() << "\t" 
                << (int) (*dfs_itr)->get_event() << "\t"
                << (*dfs_itr)->get_event_prob() << "\n";

        ++dfs_itr;
    }

    // gather outcome and latency stats (by only looking at leafs of path state 
    // tree)
    tree<Path_State *>::leaf_iterator leaf_itr      = this->path_state_tree.begin_leaf();
    tree<Path_State *>::leaf_iterator end_leaf_itr  = this->path_state_tree.end_leaf();
    while (leaf_itr != end_leaf_itr) {

        if ((modes & MODE_VERBOSE)) {

            path_state_str = (*leaf_itr)->to_string();
            std::cout << "[PATH_TREE DEPTH " 
                << (int) this->path_state_tree.depth(leaf_itr)
                << " : " << path_state_str << std::endl;

            free(path_state_str);
        }

        output_file[FILE_PATHS] << (*leaf_itr)->get_router()->get_id() << "\t"
            << (int) (*leaf_itr)->get_router()->get_height() << "\t" << (int) (*leaf_itr)->get_router()->get_width() << "\t" 
            << (int) (*leaf_itr)->get_path_status() << "\t"
            << (int) (*leaf_itr)->get_path_length() << "\t"
            << (*leaf_itr)->get_path_prob() << "\n";

        checksum += (*leaf_itr)->get_path_prob();
        ++leaf_itr;
    }

    if ((modes & MODE_VERBOSE)) 
        std::cout << "\n[CHECKSUM : " << checksum << std::endl;

    return 0;
}
