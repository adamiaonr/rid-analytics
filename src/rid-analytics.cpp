#include <math.h>
#include <cfloat>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>  // reading .nw files
#include <map>
#include <list>
#include <bitset>

#include "graph.h"
#include "rid-analytics.h"
#include "pugixml.hpp"

RID_Analytics::RID_Analytics(
    std::string nw_filename,
    uint8_t request_size,
    uint16_t bf_size,
    std::string origin_server_id,
    int mode,
    std::string output_dir,
    std::string output_label) {

    // initialize the scenario parameters set at runtime
    this->request_size = request_size;
    this->bf_size = bf_size;
    this->mode = mode;
    std::cout << "RID_Analytics::RID_Analytics() [INFO] : mode = " 
        << std::bitset<(sizeof(int) * 8)>(this->mode) << std::endl;

    // ... and build the network
    build_network(nw_filename);
    read_scn(nw_filename);

    // after building the network, save the origin server pointer
    this->origin_server = this->routers[origin_server_id];
    std::cout << "RID_Analytics::RID_Analytics() [INFO] : origin_server_id = " 
        << this->routers[origin_server_id]->get_id() << std::endl;

    // init output files
    init_output(output_label, output_dir);
}

RID_Analytics::~RID_Analytics() {

    // erase the routers in the topology. use the router map for that (see the 
    // advantage of using modern C++ features you dummy?)
    for(RID_RouterMap::iterator rit = this->routers.begin(); 
        rit != this->routers.end(); 
        ++rit) {

        delete (*rit).second;
    }
}

std::string get_unix_timestamp() {

    std::stringstream strm;
    strm << time(NULL);

    return strm.str();
}

void RID_Analytics::init_output(
    std::string label,
    std::string output_dir) {

    // we create 4 .tsv files:
    //  -# events.tsv       : probs of events (LLM, SLM, NLM, MLM) per router
    //  -# outcomes.tsv     : probs + latency + router of outcome
    //  -# forwarding.tsv   : description of fwding events, including next router probs
    //  -# ifaces.tsv       : description of fwding events, including iface probs per router

    // file names follow the convention <type>.<label>.<unix-timestamp>.tsv
    // this->output_file = std::vector<std::ofstream> (3);
    // std::string filename[3] = {"events", "outcomes", "forwarding"};
    // for (int i = 0; i < 3; i++) {
    //     filename[i] += std::string(".") + label + std::string(".") + get_unix_timestamp() + std::string(".tsv");
    //     this->output_file[i].open(output_dir + std::string("/") + filename[i]);
    // }
    this->output_file = std::vector<std::ofstream> (2);
    std::string filename[2] = {"events", "outcomes"};
    for (int i = 0; i < 2; i++) {
        filename[i] += std::string(".") + label + std::string(".") + get_unix_timestamp() + std::string(".tsv");
        this->output_file[i].open(output_dir + std::string("/") + filename[i]);
    }

    // events header
    this->output_file[0] << "router\tnlm\tmlm\tllm\tslmc\tslmi\tfpd\n";
    // outcome header
    this->output_file[1] << "router\ttype\tprob\tlatency\n";
    // // forwarding header
    // this->output_file[2] << "from\tinprob\tinto";
    // for (auto const& r : this->routers)
    //     this->output_file[2] << "\t" << r.second->get_id();
    // this->output_file[2] << "\n";
    
    // ifaces header
    // determine router w/ largest iface_num
    uint8_t iface_num = 0;
    for (auto const& r : this->routers) {
        iface_num = std::max(iface_num, r.second->get_iface_num());
    }
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

            std::vector<uint8_t> * tree_bitmask = next_router->get_tree_bitmask(k);
            // std::cout << "RID_Analytics::get_origin_distance() : [INFO] "
            //     << "\n\ttree_bitmask = " << (int) tree_bitmask[bitmask_byte] 
            //     << "\n\tbitmask_value = " << (int) bitmask_value << std::endl;

            if ((bitmask_value & (*tree_bitmask)[bitmask_byte]) == bitmask_value) {

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

    this->f_r_distribution = std::vector<__float080>(this->request_size, 0.0);
    for (pugi::xml_node f_r_dist = topology.child("f_r_dist"); 
        f_r_dist; 
        f_r_dist = f_r_dist.next_sibling("f_r_dist")) {

        int i = (f_r_dist.attribute("diff").as_int() - 1);
        this->f_r_distribution[i] = (__float080) f_r_dist.text().as_double();

        std::cout << "RID_Analytics::read_scn() : [INFO] |F\\R|[" 
            << i << "] = " << this->f_r_distribution[i] << std::endl;
    }

    // tp sizes per router
    for (pugi::xml_node router = topology.child("router"); 
        router; 
        router = router.next_sibling("router")) {

        std::string router_id = std::string(router.attribute("id").value());

        // int * _tp_sizes = (int *) calloc(routers[router_id]->get_iface_num(), sizeof(int));
        std::vector<uint8_t> _tp_sizes(routers[router_id]->get_iface_num(), 0);
        for (pugi::xml_node tp = router.child("tp"); 
            tp; 
            tp = tp.next_sibling("tp")) {

            _tp_sizes[tp.attribute("iface").as_uint()] = tp.text().as_uint();
        }

        // we keep tp sizes in a map router_id <-> tp sizes (per iface)
        this->tp_sizes[router_id] = _tp_sizes;
        this->routers[router_id]->set_tp_sizes(&(_tp_sizes));
        this->routers[router_id]->set_f_r_distr(&(this->f_r_distribution));
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
        std::cout << "RID_Analytics::build_network() : [INFO] installing router " << router_id << std::endl;

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

        std::cout << "RID_Analytics::build_network() : [INFO] iface_num = " << (int) iface_num << std::endl;

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
                this->request_size,         // max. request & forwarding entry size
                this->bf_size,              // BF size (of both requests and entries, in bit)
                this->mode);                // modes for multiple matches and error resolution

        } else {

            new_router = 
                new RID_Router(
                    router_id,
                    table_size,
                    (uint8_t) iface_num,
                    this->request_size,
                    this->bf_size,
                    this->mode);

            // add the router to the router map
            this->routers[router_id] = new_router;

            std::cout << "RID_Analytics::build_network() : [INFO] added new router : "
                << router_id << std::endl;
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

            // ***
            // *** extract the tree bitmask : a bitmask which indicates which
            // *** sources are reachable from this iface 
            // ***
            std::map<int, std::vector<uint8_t> > tree_bitmasks;

            for (pugi::xml_node tree_bitmask_block = link.child("tree_bitmask"); 
                tree_bitmask_block; 
                tree_bitmask_block = tree_bitmask_block.next_sibling("tree_bitmask")) {

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

                unsigned int tb_size = tree_bitmask_block.attribute("size").as_uint();
                tree_bitmasks[tb_size] = tree_bitmask;

                std::cout << "RID_Analytics::build_network() : [INFO] tree_bitmask["
                    << router_id << "][" << link.attribute("local").as_int() << "] : ";
                for (unsigned int i = 0; i < tree_bitmasks[tb_size].size(); i++)
                    std::cout << "[" << (int) tree_bitmasks[tb_size][i] << "] ";
                std::cout << std::endl;
            }
            
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
                &tree_bitmasks,
                adjacent_router, 
                link.attribute("remote").as_int());
        }
    }

    return 0;
}

void RID_Analytics::add_outcome(
    RID_Router * router,
    int outcome, 
    __float080 prob,
    __float080 latency) {

    // header : router\ttype\tprob\tlatency\n
    std::string code = "";
    switch (outcome) {

        case OUTCOME_CORRECT_DELIVERY:
            code = "cd";
            break;

        case OUTCOME_INCORRECT_DELIVERY:
            code = "id";
            break;

        case OUTCOME_SOFT_FALLBACK_FWD:
            code = "sff";
            break;

        case OUTCOME_SOFT_FALLBACK_AID_FWD:
            code = "sffid";
            break;

        case OUTCOME_HARD_FALLBACK_DELIVERY:
            code = "hfd";
            break;

        case OUTCOME_FEEDBACK_DELIVERY:
            code = "fed";
            break;

        case OUTCOME_PACKET_DROP:
            code = "dr";
            break;

        default:
            code = "uk";
            break;
    }

    // add line to outcome file
    this->output_file[1] << router->get_id() << "\t" 
        << code << "\t" << prob << "\t" << latency << "\n";
}

bool on_path_to_origin(RID_Router * router, int i, int os_id) {

    // check if bit is on 
    unsigned int byte_ix = (unsigned int) (os_id / 8);
    std::vector<uint8_t> * bitmask = router->get_tree_bitmask(i);

    // std::cout << "RID_Analytics : on_path_to_origin() : is [" 
    //     << router->get_id() << "][" << i << "] on the path to " << os_id << "?"
    //     << "\n\t [byte] : " << byte_ix
    //     << "\n\t [bit] : " << std::bitset<(sizeof(uint8_t) * 8)>(1 << (os_id % 8))
    //     << "\n\t [ bitmask[byte] ] : " << std::bitset<(sizeof(uint8_t) * 8)>((*bitmask)[byte_ix])
    //     << "\n\t [result] : " << (((1 << (os_id % 8)) & (*bitmask)[byte_ix]) > 0 ? "YES" : "NO") << std::endl;

    if (byte_ix < (*bitmask).size())
        return (((1 << (os_id % 8)) & (*bitmask)[byte_ix]) > 0);

    return false;
}

void RID_Analytics::add_events(
    RID_Router * router,
    std::vector<__float080> events) {

    // header: router\tnlm\tllm\tmlm\tslm\n
    this->output_file[0] << router->get_id();
    for (unsigned int e = 0; e < events.size(); e++) {

        // if(((e == EVENT_SLMC) || (e == EVENT_SLMI) || (e == EVENT_LLM)) && events[e] > 1.1)
        //     std::cout << "RID_Analytics::add_events() : [ERROR] E[" << ((e == EVENT_LLM) ? "LLM" : "SLMx") << "] = " 
        //         << events[e] << " (> 1.0) @R" << router->get_id() << std::endl;

        this->output_file[0] << "\t" << events[e];
    }
    this->output_file[0] << "\n";
}

void RID_Analytics::add_fwd(
    RID_Router * router,
    int in_iface,
    Path_State * state,
    std::vector<RID_Analytics::fwd_decision> fwd_decisions) {

    // forwarding file
    // forwarding header : 'from\tinprob\tinto\t<router-0>\t<router-1>\t....\t<router-n>\n'
    this->output_file[2] << router->get_id() 
        << "\t" << state->get_in_iface_prob() 
        << "\t" << router->get_next_hop(in_iface).router->get_id();

    // add routers attached to ifaces of router
    std::map<std::string, __float080> neighbor_probs;
    for (unsigned int i = 0; i < fwd_decisions.size(); i++) {

        std::string neighbor_id = fwd_decisions[i].next_router->get_id();
        neighbor_probs[neighbor_id] = fwd_decisions[i].state->get_in_iface_prob();
    }

    // add line to forwarding file
    for (auto const& r : this->routers) {
        if (neighbor_probs.count(r.second->get_id()) > 0)
            this->output_file[2] << "\t" << neighbor_probs[r.second->get_id()];
        else
            this->output_file[2] << "\t" << 0.0;
    }
    this->output_file[2] << "\n";
}

int RID_Analytics::handle_llm(
    RID_Router * router,
    std::vector<std::vector<__float080> > & iface_probs,
    std::vector<__float080> & events,
    Path_State * prev_state) {

    // LLM : 'local link match' i.e. a delivery to a cache
    std::cout << "RID_Analytics::handle_llm() : analyzing LLM event" 
        << std::endl;

    // update the events array
    events[EVENT_LLM] = iface_probs[0][2];

    // don't continue if prob of llm is 0.0
    if (iface_probs[0][2] == 0.0) return 0;
    
    // assess correctness of local delivery
    if (this->tp_sizes[router->get_id()][0] > 0) {

        // add record of correct delivery
        add_outcome(router, OUTCOME_CORRECT_DELIVERY, iface_probs[0][2], prev_state->get_length());

    } else {

        // add record of incorrect delivery
        add_outcome(router, OUTCOME_INCORRECT_DELIVERY, iface_probs[0][2], prev_state->get_length());

        // handling incorrect deliveries if resolution is on
        if (this->mode & 0x04) {

            if (this->mode == MMH_FLOOD) {

                add_outcome(
                    router, 
                    OUTCOME_FEEDBACK_DELIVERY,
                    iface_probs[0][2], 
                    prev_state->get_length() + get_origin_distance(this->start_router));

            } 
            // else if ((this->mode & 0x03) == HARD_FALLBACK) {

            //     add_outcome(
            //         router, 
            //         OUTCOME_HARD_FALLBACK_DELIVERY, 
            //         iface_probs[0][0],
            //         prev_state->get_length() + get_origin_distance(router));

            // }
        }
    }

    return 0;
}

int RID_Analytics::handle_ttl_drop(
    RID_Router * router,
    __float080 event_num,
    Path_State * prev_state) {

    std::cout << "RID_Analytics::handle_ttl_drop() : analyzing TTL drop event" 
        << std::endl;

    add_outcome(
        router, 
        OUTCOME_TTL_DROP, 
        event_num, 
        prev_state->get_length());

    // FIXME: special event to end simulation
    //  - it's not solved by any kind of resolution
    return 0;
}

int RID_Analytics::add_fwd_decision(
    RID_Router * router,
    Path_State * prev_state,
    std::vector<std::vector<__float080> > out_fptree_probs,
    int iface,
    __float080 path_prob,
    int latency_incr,
    std::vector<RID_Analytics::fwd_decision> & fwd_decisions) {

    Path_State * next_state = new Path_State(router, this->request_size);

    next_state->set_length(prev_state->get_length() + latency_incr);
    next_state->set_in_fptree_prob(&out_fptree_probs[iface], this->request_size);
    next_state->set_path_prob(path_prob);
    next_state->set_in_iface_prob(path_prob);
    next_state->set_ttl(prev_state->get_ttl() - 1);

    RID_Router::nw_address next_hop = router->get_next_hop(iface);
    next_state->set_tree_bitmasks(router->get_tree_bitmasks(iface));

    std::cout << "RID_Analytics::add_fwd_decision() : [INFO] saving fwd decision :"
        << "\n\tout[" << router->get_id() << "][" << iface << "] -> in[" << next_hop.router->get_id() << "][" << (int) next_hop.iface << "]"
        << "\n\tpath prob = " << path_prob
        << "\n\tpath length = " << prev_state->get_length() + latency_incr
        << "\n\tttl = " << (int) next_state->get_ttl() << std::endl;

    RID_Analytics::fwd_decision fd = {router, next_hop.router, iface, next_hop.iface, next_state};
    fwd_decisions.push_back(fd);

    return 0;
}

int RID_Analytics::make_fwd_decision(
    RID_Router * router,
    Path_State * prev_state,
    uint8_t in_iface,
    std::vector<std::vector<__float080> > iface_probs,
    std::vector<std::vector<__float080> > out_fptree_probs,
    std::vector<RID_Analytics::fwd_decision> & fwd_decisions) {

    // IMPORTANT: here's what the 3 values in iface_probs[i][] carry
    //  - iface_probs[i][0] : P(L_i >  L_{~i}), i.e. the probability of an 'exclusive'
    //                        match at iface i.
    //                        this should be used for Single Match (SM) events
    //  - iface_probs[i][1] : P(L_i >= L_{~i}), Multiple Match (MM) situations which 
    //                        *DO NOT* involve the choice of iface 0.
    //                        this is useful when evaluating the probability of 
    //                        'fallback' events, in which we want to distinguish 
    //                        MM situations which involve / not involve local 
    //                        matches.
    //  - iface_probs[i][2] : P(L_i >= L_{~i}), multiple match situations which 
    //                        consider the involvement of iface 0.
    //                        this should be the maximum value of all iface[i][].

    // array to save avg. nr of events
    std::vector<__float080> events(EVENT_NUM, 0.0);

    // handle LLM events (iface 0)
    handle_llm(router, iface_probs, events, prev_state);

    // FIXME: special vars for soft fallback handling
    __float080 max_prob[3] = {0.0, 0.0, 0.0};
    int origin_iface = 0;
    // handle other events (iface > 0)
    for (int i = 1; i < router->get_iface_num(); i++) {

        if (on_path_to_origin(router, i, std::stoi(this->origin_server->get_id())))
            origin_iface = i;

        // check if the packet's ttl goes below 0 : if it does, end the path here
        // without adding a forwarding decision
        if ((prev_state->get_ttl() - 1) < 0) { 
            if (iface_probs[i][2] > 0.0)
                handle_ttl_drop(router, iface_probs[i][2], prev_state);
            continue;
        }

        // update the path probability (P(R_n) in the paper)
        __float080 path_prob = 0.0;
        if (this->mode == MMH_FLOOD) {

            // MMH_FLOOD uses P(L_i >= L_{~i}), w/ local matches
            path_prob = iface_probs[i][2];
            // SM event only considers P(L_i > L_{~i})

            if (this->tp_sizes[router->get_id()][i] > 0) {
                events[EVENT_SLMC] += iface_probs[i][0];
            } else {
                events[EVENT_SLMI] += iface_probs[i][0];
            }

            // MM event refers to P(L_i == L_{~i}), as such we do:
            //  P(L_i >= L_{~i}) - P(L_i > L_{~i})
            // FIXME : this will be used to calculate the avg. nr. of links 
            // used in 'flood' mode
            events[EVENT_MLM] += iface_probs[i][2] - iface_probs[i][0];

        } else if (((this->mode & 0x03) == HARD_FALLBACK)
            || ((this->mode & 0x03) == SOFT_FALLBACK)) {

            // as an approximation for the MM event probability, 
            // we use the iface which outputs the max. value of iface_probs[i][1]
            // FIXME: how do we explain this?
            if (max_prob[1] < iface_probs[i][2]) {

                max_prob[0] = iface_probs[i][0];
                max_prob[1] = iface_probs[i][1];

                // FIXME: workaround to get correct value of max_prob at start router
                if ((*router->get_blocked_ifaces())[0])
                    max_prob[2] = iface_probs[i][1];
                else
                    max_prob[2] = iface_probs[i][2];
            }

            // otherwise, use iface_probs[i][0]
            path_prob = iface_probs[i][0];
            if (this->tp_sizes[router->get_id()][i] > 0)
                events[EVENT_SLMC] += iface_probs[i][0];
            else
                events[EVENT_SLMI] += iface_probs[i][0];
        }

        // if path prob = 0.0, skip the addition of a fwd decision
        if (path_prob == 0.0) continue;

        // create a new path state
        Path_State * next_state = new Path_State(router, this->request_size);
        next_state->set_length(prev_state->get_length() + 1);
        // now, we deal with probabilities:
        //
        //  -# in_fptree_probs : the prob of having the packet 
        //     bound to a prefix tree of size s (no TP info)
        //
        //  -# in_iface_probs : the prob of having a packet 
        //     flow through iface i (takes TPs into account)
        //
        next_state->set_in_fptree_prob(&out_fptree_probs[i], this->request_size);
        next_state->set_path_prob(path_prob);
        next_state->set_in_iface_prob(path_prob);
        next_state->set_ttl(prev_state->get_ttl() - 1);

        // finally, we continue with the request path get the next hop router by iface
        RID_Router::nw_address next_hop = router->get_next_hop(i);
        // set the path's tree bitmask (by extracting the bitmask from the router's iface)
        next_state->set_tree_bitmasks(router->get_tree_bitmasks(i));

        // this can still happen (or can it?)
        // FIXME: this subtle condition is the one that actually terminates the run. 
        // FIXME: this could be the source of trouble.
        if (next_hop.router == router || i == in_iface) {
            std::cout << "RID_Analytics::make_fwd_decision() : [INFO] "
                << "found magic stopping condition" << std::endl;

            delete next_state;
            continue;
        }

        // add fwd_decision to fwd decision vector
        std::cout << "RID_Analytics::make_fwd_decision() : [INFO] saving fwd decision :"
            << "\n\tout[" << router->get_id() << "][" << i << "] -> in[" << next_hop.router->get_id() << "][" << (int) next_hop.iface << "]"
            << "\n\tpath prob = " << path_prob
            << "\n\tttl = " << (int) next_state->get_ttl() << std::endl;

        RID_Analytics::fwd_decision fd = {router, next_hop.router, i, next_hop.iface, next_state};
        fwd_decisions.push_back(fd);
    }

    // FIXME: special handling of SOFT_FALLBACK case
    if (this->mode == MMH_FLOOD) {

        events[EVENT_FPD] = std::max(
            (__float080) 0.0, 
            (prev_state->get_in_iface_prob() - events[EVENT_SLMC] - events[EVENT_SLMI] - iface_probs[0][0]));

    } else if ((this->mode & 0x03) == SOFT_FALLBACK) {

        bool updated = false;
        for (unsigned int k = 0; k < fwd_decisions.size(); k++) {

            if (fwd_decisions[k].out_iface == origin_iface) {

                // update fwd_decision on fwd decision vector
                fwd_decisions[k].state->set_path_prob(max_prob[1]);
                fwd_decisions[k].state->set_in_iface_prob(max_prob[1]);
                std::cout << "RID_Analytics::make_fwd_decision() : [INFO] updating fwd decision :"
                    << "\n\tout[" << router->get_id() << "][" << origin_iface << "] -> in[" << fwd_decisions[k].next_router->get_id() << "][" << (int) fwd_decisions[k].in_iface << "]"
                    << "\n\tpath prob = " << max_prob[1]
                    << "\n\tpath length = " << fwd_decisions[k].state->get_length()
                    << "\n\tttl = " << (int) fwd_decisions[k].state->get_ttl() << std::endl;

                // mark fwd decision over origin iface as 'updated'
                updated = true;

                // add 'soft fallback fwd' outcome w/ probability P(L_i == L_{~i})
                if ((max_prob[1] - max_prob[0]) > 0.0)
                    add_outcome(
                        router, 
                        OUTCOME_SOFT_FALLBACK_FWD, 
                        (max_prob[1] - max_prob[0]), 
                        prev_state->get_length());

                events[EVENT_MLM] += (max_prob[2] - max_prob[0]);
            }
        }

        // in case a fwd decision isn't updated above, add a new fwd decision
        if ((!updated) && (max_prob[1] > 0.0)) {
            add_fwd_decision(router, prev_state, out_fptree_probs, origin_iface, max_prob[1], 1, fwd_decisions);
        }

        // handle the case of having the request incorrectly delivered and 
        // recovered w/ a soft fallback (w/ +2 latency penalty)
        if ((!(this->tp_sizes[router->get_id()][0] > 0)) && ((max_prob[2] - max_prob[1]) > 0.0) && (iface_probs[0][2] > 0.0)) {
            add_fwd_decision(router, prev_state, out_fptree_probs, origin_iface, max_prob[2] - max_prob[1], 2, fwd_decisions);
        }

    } else if ((this->mode & 0x03) == HARD_FALLBACK) {

        __float080 mm_prob = std::max(
            (__float080) 0.0, 
            (prev_state->get_in_iface_prob() - events[EVENT_SLMC] - events[EVENT_SLMI] - iface_probs[0][2]));
        
        events[EVENT_FPD] = mm_prob;

        if (mm_prob > 0.0) {

            // 1) add correct delivery 
            add_outcome(
                router, 
                OUTCOME_CORRECT_DELIVERY, 
                (mm_prob), 
                prev_state->get_length() + get_origin_distance(router));

            // 2) register hard fallback activation
            add_outcome(
                router, 
                OUTCOME_HARD_FALLBACK_DELIVERY, 
                (mm_prob), 
                get_origin_distance(router));

            events[EVENT_MLM] += mm_prob;
        }

        if ((iface_probs[0][2] > 0.0) && (router->get_id() != origin_server->get_id()) && (!(this->tp_sizes[router->get_id()][0] > 0))) {

            events[EVENT_SLMI] += iface_probs[0][2];

            // 1) add correct delivery 
            add_outcome(
                router, 
                OUTCOME_CORRECT_DELIVERY, 
                iface_probs[0][2], 
                prev_state->get_length() + 1 + get_origin_distance(router));

            // 2) register hard fallback activation
            add_outcome(
                router, 
                OUTCOME_HARD_FALLBACK_DELIVERY, 
                iface_probs[0][2], 
                get_origin_distance(router));
        }
    }

    // // add fwd record to results
    // add_fwd(router, in_iface, prev_state, fwd_decisions);
    //  - add events probs to results
    add_events(router, events);

    return 0;
}

int RID_Analytics::run_rec(
    RID_Router * router,
    uint8_t in_iface,
    Path_State * prev_state) {

    std::cout << "\nRID_Analytics::run_rec() : forward operation by router : " 
        << router->get_id() << std::endl;

    uint8_t iface_num = router->get_iface_num();
    __float080 in_prob = prev_state->get_in_iface_prob();
    std::vector<__float080> * in_fptree_prob = prev_state->get_in_fptree_prob();
    // tree bitmask of the path tells us which sources are accessible over 
    // the followed path.
    std::map<int, std::vector<uint8_t> > * tree_bitmasks = prev_state->get_tree_bitmasks();

    // run the forward() operation on this router & collect: 
    //  - iface probs
    //  - event probs
    //  - out fp tree probs 
    std::vector<std::vector<__float080> > iface_probs(iface_num, std::vector<__float080> (3, 0.0));
    std::vector<std::vector<__float080> > out_fptree_probs(iface_num, std::vector<__float080> (this->request_size + 1, 0.0));
    router->forward(
        in_iface,           // INPUTS
        tree_bitmasks,
        in_prob,
        in_fptree_prob,
        iface_probs,        // OUTPUTS
        out_fptree_probs);

    // make fwd decisions
    std::vector<RID_Analytics::fwd_decision> fwd_decisions;
    make_fwd_decision(router, prev_state, in_iface, iface_probs, out_fptree_probs, fwd_decisions);
    // recursively call run_rec() to resume fwd decisions
    for (unsigned int k = 0; k < fwd_decisions.size(); k++) {

        std::cout << "RID_Analytics::run_rec() : [INFO] on router[" << router->get_id() 
            << "], forwarding to router[" << fwd_decisions[k].next_router->get_id() 
            << "], ttl = " << (int) fwd_decisions[k].state->get_ttl() << std::endl;

        run_rec(fwd_decisions[k].next_router, fwd_decisions[k].in_iface, fwd_decisions[k].state);
    }

    delete prev_state;

    return 0;
}

void RID_Analytics::run(
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

    // FIXME: not sure if it is ok to set the rid_router arg on Path_State() as NULL
    Path_State * init_state = new Path_State(NULL, request_size);
    init_state->set_length(0);

    // set the initial ingress probabilities: ptree_size = 0 gets 1.0 
    // probability, i.e. initially the request isn't bound to any prefix tree
    init_state->set_in_fptree_prob(0, 1.0);
    init_state->set_in_iface_prob(1.0);

    // all other prefix tree sizes get 0.0 probability
    for (uint8_t p = 1; p <= this->request_size; p++)
        init_state->set_in_fptree_prob(p, 0.0);

    // we first supply an empty tree bitmask for the path
    std::map<int, std::vector<uint8_t> > tree_bitmasks =
        (*this->start_router->get_tree_bitmasks(0));
    // fill the path tree bitmask w/ 0s
    std::map<int, std::vector<uint8_t> >::iterator itr;
    for (itr = tree_bitmasks.begin(); itr != tree_bitmasks.end(); itr++)
        std::fill(itr->second.begin(), itr->second.end(), 0x00);

    // set the path tree bitmask on path state
    init_state->set_tree_bitmasks(&tree_bitmasks);
    // set initial ttl to dist. from starting router to the origin server
    init_state->set_ttl(4);

    // we always start at the ingress iface, so that we don't have a local 
    // match at the initial router (not sure if this can work directly like this)
    run_rec(this->start_router, 0, init_state);
}
