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
#include "pugixml.hpp"

RID_Analytics::RID_Analytics(
    std::string nw_filename,
    uint8_t request_size,
    uint16_t bf_size,
    std::string origin_server_id,
    int mm_mode,
    int resolution_mode) {

    // initialize the scenario parameters set at runtime
    this->request_size = request_size;
    this->f_max = request_size;
    this->bf_size = bf_size;
    // modes:
    //  -# mm_mode : modes for handling multiple matches in routers. this 
    //               influences efficiency. there are 3 diff. modes:
    //                  -# MMH_FLOOD    : forward over all matching ifaces
    //                  -# MMH_RANDOM   : forward over 1 random iface (not used)
    //                  -# MMH_FALLBACK : forward using the fallback address
    //
    //
    //  -# resolution_mode : specifies if wrong deliveries should be fixed 
    //                       or not. it influences the probability of (in)correct
    //                       deliveries and latency. it can take 2 values:
    //                          -# RESOLUTION_OFF : no fixes if a packet is 
    //                             delivered incorrectly.
    //                          -# RESOLUTION_ON : packet is sent to the origin 
    //                             server if a wrong delivery occurs.  
    //                             
    this->mm_mode = mm_mode;
    this->resolution_mode = resolution_mode;

    this->f_r_distribution = NULL;

    // ... and build the network
    this->build_network(nw_filename);
    this->read_scn(nw_filename);

    // after building the network, save the origin server pointer
    this->origin_server = this->routers[origin_server_id];
    std::cout << "RID_Analytics::RID_Analytics() [INFO] : origin_server_id = " 
        << this->routers[origin_server_id]->get_id() << std::endl;

    // // FIXME : hack to make simulations shorter
    // this->hit = false;
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

int RID_Analytics::on_path_to_origin(
    RID_Router * router, 
    int ingress_iface, 
    int egress_iface = -1) {

    uint8_t bitmask_value = (1 << (std::stoi(this->origin_server->get_id()) % 8));
    uint8_t bitmask_byte = (std::stoi(this->origin_server->get_id()) / 8);
    uint8_t tree_bitmask_value = 0;

    if ((egress_iface == -1)) {

        for (uint8_t k = 1; k < router->get_iface_num(); k++) {

            if (k == ingress_iface)
                continue;

            tree_bitmask_value = router->get_tree_bitmask(k)[bitmask_byte];
            if ((bitmask_value & tree_bitmask_value) == bitmask_value)
                return 1;
        }

        return 0;
    }

    tree_bitmask_value = router->get_tree_bitmask(egress_iface)[bitmask_byte];
    if ((bitmask_value & tree_bitmask_value) == bitmask_value)
        return 1;

    return 0;
}

int RID_Analytics::get_origin_distance(RID_Router * from_router) {

    // check if we're there yet
    if (from_router->get_id() == this->origin_server->get_id())
        return 0;

    uint8_t bitmask_value = (1 << (std::stoi(this->origin_server->get_id()) % 8));
    uint8_t bitmask_byte = (std::stoi(this->origin_server->get_id()) / 8);

    // std::cout << "RID_Analytics::get_origin_distance() : [INFO] "
    //     << "\n\torigin_server_id = " << this->origin_server->get_id() 
    //     << "\n\tbitmask_value = " << (int) bitmask_value 
    //     << "\n\tbitmask_byte = " << (int) bitmask_byte << std::endl;

    int origin_distance = 0;
    RID_Router * next_router = from_router;
    while (next_router->get_id() != this->origin_server->get_id()) {

        // std::cout << "RID_Analytics::get_origin_distance() : [INFO] "
        //     << "\n\tnext router id = " << next_router->get_id() << std::endl;

        // get next router in path
        uint8_t iface_num = next_router->get_iface_num();
        for (uint8_t k = 0; k < iface_num; k++) {

            uint8_t * tree_bitmask = next_router->get_tree_bitmask(k);
            // std::cout << "RID_Analytics::get_origin_distance() : [INFO] "
            //     << "\n\ttree_bitmask = " << (int) tree_bitmask[bitmask_byte] 
            //     << "\n\tbitmask_value = " << (int) bitmask_value << std::endl;

            if ((bitmask_value & tree_bitmask[bitmask_byte]) == bitmask_value) {

                // save the next router in the path in next_router
                next_router = next_router->get_next_hop(k).router;
                // increment the origin distance
                origin_distance++;
                // jump off the cycle
                k = iface_num;
            }
        }
    }

    return origin_distance;
}

int RID_Analytics::read_scn(std::string nw_filename) {

    // extract the TP size map and |F\R| distr. from nw_filename
    // read .scn file, structured as an xml file
    pugi::xml_document nw_doc;
    if (!(nw_doc.load_file(nw_filename.c_str()))) {

        std::cerr << "RID_Analytics::get_network() : [ERROR] .scn file (" 
            << nw_filename << ") not found." << std::endl;

        return -1;
    }

    pugi::xml_node topology = nw_doc.child("topology");

    // tp sizes per router
    for (pugi::xml_node router = topology.child("router"); 
        router; 
        router = router.next_sibling("router")) {

        std::string router_id = std::string(router.attribute("id").value());

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
    for (pugi::xml_node f_r_dist = topology.child("f_r_dist"); 
        f_r_dist; 
        f_r_dist = f_r_dist.next_sibling("f_r_dist")) {

        int i = (f_r_dist.attribute("diff").as_int() - 1);
        this->f_r_distribution[i] = (__float080) f_r_dist.text().as_double();

        std::cout << "RID_Analytics::read_scn() : [INFO] |F\\R|[" 
            << i << "] = " << this->f_r_distribution[i] << std::endl;
    }

    return 0;
}

int RID_Analytics::build_network(std::string nw_filename) {

    // read .scn file, structured as an xml file
    pugi::xml_document nw_doc;
    if (!(nw_doc.load_file(nw_filename.c_str()))) {

        std::cerr << "RID_Analytics::get_network() : [ERROR] .scn file (" 
            << nw_filename << ") not found." << std::endl;

        return -1;
    }

    pugi::xml_node topology = nw_doc.child("topology");

    // ***
    // *** extract the fwd table size (a unique table size for the complete topology)
    // ***
    std::cout << "RID_Analytics::build_network() : [INFO] table size: ";
    pugi::xml_node ts = topology.child("fwd_table_size");
    uint64_t table_size = ts.text().as_uint();
    std::cout << "\n\tSIZE = " << ts.text().as_uint() << " (" << table_size << ")" << std::endl;

    // ***
    // *** extract the ttl
    // ***
    std::cout << "RID_Analytics::build_network() : [INFO] ttl: ";
    pugi::xml_node ttl_block = topology.child("ttl");
    this->ttl = ttl_block.text().as_int();
    std::cout << "\n\tTTL = " << ttl_block.text().as_int() << " (" << ttl << ")" << std::endl;

    // ***
    // *** start building the topology by reading the info for each router
    // ***
    for (pugi::xml_node router = topology.child("router"); 
        router; 
        router = router.next_sibling("router")) {

        std::string router_id = std::string(router.attribute("id").value());

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
                router_id,
                table_size,                 // nr. of table entries for router
                iface_num,                  // nr. of ifaces (i.e. nr. of adjacent routers)
                this->f_max,                // max. request & forwarding entry size
                this->bf_size,              // BF size (of both requests and entries, in bit)
                this->mm_mode);             // multiple match mode

        } else {

            new_router = 
                new RID_Router(
                    router_id,
                    table_size,
                    iface_num,
                    this->f_max,
                    this->bf_size,
                    this->mm_mode);

            // add the router to the router map
            this->routers[router_id] = new_router;

            std::cout << "RID_Analytics::build_network() : [INFO] added new router : "
                << router_id << std::endl;
        }

        // if (tier == this->access_tree_height && tier_index > 0)
        //     new_router->set_leaf();

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

            // ***
            // *** extract the tree bitmask : a bitmask which indicates which
            // *** sources are reachable from this iface 
            // ***
            pugi::xml_node tree_bitmask_block = link.child("tree_bitmask");
            std::string tree_bitmask_str = std::string(tree_bitmask_block.text().as_string());

            std::vector<uint8_t> tree_bitmask;
            for (unsigned int i = 0; i < tree_bitmask_str.size(); ) {

                uint8_t higher_bits = (uint8_t) ((tree_bitmask_str[i] >= 'a') ? (tree_bitmask_str[i] - 'a' + 10) : (tree_bitmask_str[i] - '0'));
                // std::cout << "RID_Analytics::build_network() : [INFO] tree_bitmask["
                //     << router_id << "][" << link.attribute("local").as_int() << "] : "
                //     << "\n\tstr (H4): " << tree_bitmask_str[i]
                //     << "\n\tint (H4): " << (int) higher_bits << std::endl;

                i++;

                uint8_t lower_bits = (uint8_t) ((tree_bitmask_str[i] >= 'a') ? (tree_bitmask_str[i] - 'a' + 10) : (tree_bitmask_str[i] - '0'));
                // std::cout << "RID_Analytics::build_network() : [INFO] tree_bitmask["
                //     << router_id << "][" << link.attribute("local").as_int() << "] : "
                //     << "\n\tstr (L4): " << tree_bitmask_str[i]
                //     << "\n\tint (L4): " << (int) lower_bits
                //     << "\n\tint (H4L4): " << (int) (((higher_bits << 4) & 0xF0) | (lower_bits & 0x0F)) << std::endl;

                i++;

                tree_bitmask.push_back(((higher_bits << 4) & 0xF0) | (lower_bits & 0x0F));
            }

            std::cout << "RID_Analytics::build_network() : [INFO] tree_bitmask["
                << router_id << "][" << link.attribute("local").as_int() << "] : ";
            for (unsigned int i = 0; i < tree_bitmask.size(); i++)
                std::cout << "[" << (int) tree_bitmask[i] << "] ";
            std::cout << std::endl;
            
            std::string adjacent_router_id = std::string(link.attribute("rrouter").value());

            // if adjacent_router is already in the router map, fetch its 
            // pointer. otherwise, create a new RID_Router object.
            RID_Router * adjacent_router = NULL;
            RID_RouterMap::iterator rtr_map_it;       
            if ((rtr_map_it = this->routers.find(adjacent_router_id)) != this->routers.end()) {

                adjacent_router = (*rtr_map_it).second;

            } else {

                adjacent_router = new RID_Router();
                // add adj. router to the router map
                this->routers[adjacent_router_id] = adjacent_router;

                std::cout << "RID_Analytics::build_network() : [INFO] added new router : "
                    << adjacent_router_id << std::endl;
            }

            new_router->add_fwd_table_entry(
                link.attribute("local").as_int(),
                fwd_dist[link.attribute("local").as_uint()],
                &size_dist,
                &tree_bitmask,   // NEW
                adjacent_router, 
                link.attribute("remote").as_int());
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

    __float080 path_prob = 0.0, fallback_carry_prob = 0.0;
    __float080 ingress_prob = (*prev_path_state_itr)->get_ingress_iface_prob();
    __float080 * ingress_ptree_prob = (*prev_path_state_itr)->get_ingress_ptree_prob();
    // tree bitmask of the path tells us which sources are accessible over 
    // the followed path.
    // FIXME : we assume the tree bitmask size is the same at every iface on 
    // the topology 
    uint8_t * tree_bitmask = (*prev_path_state_itr)->get_tree_bitmask();

    // next hop information
    RID_Router::nw_address next_hop;
    tree<Path_State *>::iterator path_state_itr;
    std::vector<RID_Analytics::run_record> single_link_match_states;

    // run the forward() operation on this router : this will determine the 
    // following info:
    //  1) interface event probabilities
    //  2) egress size probabilities
    //
    // we then use this info to determine the path state
    std::cout << "\nRID_Analytics::run_rec() : FORWARD() BY ROUTER : " 
        << router->get_id() << std::endl;

    // if (this->tp_sizes.empty())
    //     std::cout << "RID_Analytics::run_rec() : [INFO] TP sizes is empty!";
    // else {

    //     std::cout << "RID_Analytics::run_rec() : [INFO] TP sizes keys : " << std::endl;
    //     RID_TPMap::iterator itr;
    //     for (itr = this->tp_sizes.begin(); itr != this->tp_sizes.end(); itr++)
    //         std::cout << "\tTP[" << itr->first << "]" << std::endl;
    // }

    // for (int i = 0; i < routers[router->get_id()]->get_iface_num(); i++)
    //     std::cout << "\tTP[" << i << "] = " << (int) this->tp_sizes[router->get_id()][i] << std::endl;       

    router->forward(
        this->request_size,
        ingress_iface,
        this->tp_sizes[router->get_id()],
        tree_bitmask,
        ingress_prob,
        ingress_ptree_prob,
        this->f_r_distribution);

    // we now cycle through all possible interface events and set the possible 
    // path states and probabilities
    for (int event = EVENT_NLM; event < EVENT_NUM; event++) {

        // we differentiate the establishment of states by iface events
        //
        //  1) if EVENT_SLM, the path continues, and so we have 
        //     additional work : setting FP tree probabilities, 
        //     determining the intermediate state (TP or FP), calling run_rec() 
        //     again, etc.
        //
        //  2) otherwise, the path may or may not end here
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
                << "\n\tEVENT : " << event << " EVENT PROB : " << path_state->get_event_prob()
                << "\n\tPROB : " << router->get_iface_events_prob(event)
                << "\n\tPATH LATENCY : " << (*prev_path_state_itr)->get_path_length() + 1
                << std::endl;

            if (event > EVENT_NLM) {

                // why???? 
                path_state->set_eop();

                // resolving Multiple Link Match (MLM) events:
                //  1) MMH_FALLBACK : send the request over the correct iface 
                //                    towards the origin server
                //  2) MMH_FLOOD    : we skip any further processing here, and 
                //                    and handle it when processing EVENT_SLM
                if (event == EVENT_MLM) {

                    if (this->mm_mode == MMH_FALLBACK) {

                        path_state->set_path_status(OUTCOME_FALLBACK_DELIVERY);

                        if (on_path_to_origin(router, ingress_iface)) {

                            // if the router is on the path to origin, we save 
                            // the MLM event probability for later use while 
                            // processing SLM events
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
                        // this->hit = true;

                    } else {

                        path_state->set_path_status(OUTCOME_INCORRECT_DELIVERY);
                        path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

                        int penalty = 0;
                        if (this->mm_mode == MMH_FALLBACK) {

                            if (on_path_to_origin(router, ingress_iface)) {

                                // if the router is on the path to origin, we save 
                                // the LLM event probability for later use while 
                                // processing SLM events
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
                                + get_origin_distance(this->start_router);

                            if (this->resolution_mode == RESOLUTION_ON) {

                                Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                                this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                                resolution_state->set_path_status(OUTCOME_CORRECT_DELIVERY);
                                resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                                resolution_state->set_path_prob(router->get_iface_events_prob(event));
                                resolution_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1 + penalty);
                            }
                        }
                    }
                }

            } else {

                // if an NLM event, it is just dropped. some resolution 
                // method may kick in (e.g. fallbacks or feedback).
                path_state->set_path_status(OUTCOME_PACKET_DROP);
                path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

                int penalty = 0;
                if (this->mm_mode == MMH_FALLBACK) {

                    if (on_path_to_origin(router, ingress_iface)) {

                        // if the router is on the path to origin, we save 
                        // the LLM event probability for later use while 
                        // processing SLM events
                        fallback_carry_prob += router->get_iface_events_prob(event);
                        // ... and set the path length of this state to 0.0
                        // FIXME: this shouldn't be necessary, now to calculate 
                        // avg. latencies you should only look into correct 
                        // deliveries
                        path_state->set_path_status(OUTCOME_FALLBACK_DELIVERY);
                        path_state->set_path_length(0.0);

                    } else {

                        penalty = get_origin_distance(router);

                        Path_State * resolution_state = new Path_State(this->origin_server, this->request_size);
                        this->path_state_tree.append_child(prev_path_state_itr, resolution_state);

                        resolution_state->set_path_status(OUTCOME_CORRECT_DELIVERY);
                        resolution_state->set_event(EVENT_LLM, router->get_iface_events_prob(event));
                        resolution_state->set_path_prob(router->get_iface_events_prob(event));
                        resolution_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1 + penalty);
                    }

                } else {

                    penalty = (*prev_path_state_itr)->get_path_length() 
                                + get_origin_distance(this->start_router);

                    if (this->resolution_mode == RESOLUTION_ON) {

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

            // if the packet's TTL goes below 0, the path ends here. each  
            // remaining event will be translated into a NLM event.
            if (((*prev_path_state_itr)->get_ttl() - 1) < 0) {

                std::cout << "RID_Analytics::run_rec() : [INFO] TTL LIMIT REACHED" << std::endl;

                Path_State * path_state = new Path_State(router, this->request_size);
                path_state_itr = this->path_state_tree.append_child(prev_path_state_itr, path_state);

                path_state->set_path_status(OUTCOME_TTL_DROP);
                path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);
                path_state->set_event(EVENT_TTL, router->get_iface_events_prob(event));
                path_state->set_path_prob(router->get_iface_events_prob(event) + fallback_carry_prob);

                // don't do anything else
                continue;
            }

            // if EVENT_SLM, the request paths continue to be 
            // built. therefore, we cycle through all egress interfaces
            for (int iface = IFACE_LOCAL + 1; iface < iface_num; iface++) {

                // don't add a state for a blocked iface
                std::set<uint8_t> blocked_ifaces = router->get_blocked_ifaces();
                std::set<uint8_t>::iterator it = blocked_ifaces.find((uint8_t) iface);
                if (it != blocked_ifaces.end())
                    continue;

                // create a new PathState object for each egress iface
                Path_State * path_state = new Path_State(router, this->request_size);
                path_state_itr = this->path_state_tree.append_child(prev_path_state_itr, path_state);

                // increment the path length by +1 hop
                path_state->set_path_length((*prev_path_state_itr)->get_path_length() + 1);

                // now, we deal with probabilities:
                //
                //  -# ingress_ptree_probs : the prob of having the packet 
                //     bound to a prefix tree of size s (no TP info)
                //
                //  -# ingress_iface_probs : the prob of having a packet 
                //     flow through iface i (takes TPs into account)
                //
                // path_state->set_ingress_ptree_prob(router->get_egress_iface_probs(iface), this->f_max);
                path_state->set_ingress_ptree_prob(router->get_egress_ptree_prob(iface), this->f_max);

                path_prob = router->get_egress_iface_prob(iface); 
                if (on_path_to_origin(router, ingress_iface, iface) && (this->mm_mode == MMH_FALLBACK)) {

                    std::cout << "RID_Analytics::run_rec() : [INFO] router " 
                        << router->get_id() << " is on the path to the origin "
                        << ", adding " << fallback_carry_prob << " to P^E(SLM)[" 
                        << iface << "] = " << path_prob << std::endl;

                    path_prob += fallback_carry_prob;
                }

                path_state->set_path_prob(path_prob);
                path_state->set_ingress_iface_prob(path_prob);
                path_state->set_ttl((*prev_path_state_itr)->get_ttl() - 1);

                // add a record of an SLM event to the path state and its 
                // probability
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
                    << "\n\tTTL : " << (*prev_path_state_itr)->get_ttl() - 1
                    << std::endl;

                // determine the correctness of the forwarding decision
                if (this->tp_sizes[router->get_id()][iface] > 0) {

                    // if router has a TP on iface, set the even as 
                    // an intermediate TP, i.e. "we're on the right path"
                    path_state->set_path_status(STATUS_TP);
                    // if (!(this->hit))
                    //     path_state->set_ttl((*prev_path_state_itr)->get_ttl());

                } else {

                    // else, an 'unlucky' FP happened: we're on the wrong path 
                    // now (a FP match occurred, in an iface which does not 
                    // have TP entries)
                    path_state->set_path_status(STATUS_FP);
                }

                // finally, we continue with the request path
                // get the next hop router by iface
                next_hop = router->get_next_hop(iface);

                // set the path's tree bitmask (by extracting the bitmask from 
                // the router's iface)
                int _tree_bitmask_size = router->get_tree_bitmask_size(iface);
                uint8_t * _tree_bitmask = router->get_tree_bitmask(iface);
                path_state->set_tree_bitmask(_tree_bitmask, _tree_bitmask_size);

                // this can still happen (or can it?)
                // FIXME: this subtle condition is the one that actually 
                // terminates the run. 
                // FIXME: this could be the source of trouble.
                if (next_hop.router == router || iface == ingress_iface) {

                    std::cout << "RID_Analytics::run_rec() : [INFO] FOUND THE MAGIC STOPPING CONDITION" << std::endl;

                    continue;
                }

                std::cout << "RID_Analytics::run_rec() : [INFO] saving slm record :"
                    << "\n\ton router[" << router->get_id() << "]"
                    << "\n\tforwarding to router[" << next_hop.router->get_id() << "]"
                    << "\n\tttl = " << (int) path_state->get_ttl() << std::endl;


                RID_Analytics::run_record slm_record;

                slm_record.next_router = next_hop.router;
                slm_record.ingress_iface = next_hop.iface;
                slm_record.prev_path_state_itr = path_state_itr;

                single_link_match_states.push_back(slm_record);
                //run_rec(next_hop.router, next_hop.iface, path_state_itr);
            }

            std::vector<RID_Analytics::run_record>::iterator itr;
            for (itr = single_link_match_states.begin(); itr != single_link_match_states.end(); itr++) {

                std::cout << "RID_Analytics::run_rec() : [INFO] on router[" << router->get_id() 
                    << "], forwarding to router[" << (*itr).next_router->get_id() 
                    << "], ttl = " << (int) (*((*itr).prev_path_state_itr))->get_ttl() << std::endl;
                run_rec((*itr).next_router, (*itr).ingress_iface, (*itr).prev_path_state_itr);
            }
        }
    }

    return 0;
}

int RID_Analytics::run(
    std::string nw_filename,
    std::string start_router_id) {

    if (this->tp_sizes.empty())
        std::cout << "RID_Analytics::run() : [INFO] TP sizes is empty!";
    else {

        std::cout << "RID_Analytics::run() : [INFO] TP sizes keys : " << std::endl;
        RID_TPMap::iterator itr;
        for (itr = this->tp_sizes.begin(); itr != this->tp_sizes.end(); itr++)
            std::cout << "\tTP[" << itr->first << "]" << std::endl;
    }

    this->start_router = this->routers[start_router_id];

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

    // we first supply an empty tree bitmask (the size must be extracted from 
    // the first router)
    int tree_bitmask_size = this->start_router->get_tree_bitmask_size(0);
    uint8_t * tree_bitmask = (uint8_t *) calloc(tree_bitmask_size, sizeof(uint8_t));
    initial_state->set_tree_bitmask(tree_bitmask, tree_bitmask_size);
    // the initial ttl is set to the distance from starting router to the 
    // origin server
    initial_state->set_ttl(this->ttl);

    // we always start at the ingress iface, so that we don't have a local 
    // match at the initial router (not sure if this can work directly like 
    // this)
    run_rec(this->start_router, 0, path_state_tree_itr);

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
    //  -# events.tsv       : probabilities of events (LLM, SLM, NLM, MLM) per router
    //  -# outcomes.tsv     : outcomes (correct/incorrect delivery, etc.) per router
    //  -# latencies.tsv    : final path latencies (nr. of hops)
    ofstream output_file[2];  
    // file names follow the convention <type>.<label>.<unix-timestamp>.tsv
    std::string output_filename[2] = {"events", "path"};
    for (int i = 0; i < 2; i++) {
        output_filename[i] += std::string(".") + output_label + std::string(".") + get_unix_timestamp() + std::string(".tsv");
        output_file[i].open(output_dir + std::string("/") + output_filename[i]);
    }

    // add column lines to .tsv files
    output_file[FILE_EVENTS] << "AS\tEVENT\tPROB\n";
    output_file[FILE_PATHS] << "AS\tSTATUS\tLATENCY\tPROB\n";

    // gather event stats (DFS iterator traverses the whole path state tree)
    tree<Path_State *>::iterator dfs_itr = this->path_state_tree.begin();
    while(dfs_itr != this->path_state_tree.end()) {

        if ((*dfs_itr)->get_router())
            output_file[FILE_EVENTS] << (*dfs_itr)->get_router()->get_id() << "\t"
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
