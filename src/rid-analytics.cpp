#include <math.h>
#include <cfloat>
#include <algorithm>
#include <string>

#include "graph.h"
#include "rid-analytics.h"

RID_Analytics::RID_Analytics(
    uint8_t access_tree_num,
    uint8_t access_tree_height,
    uint8_t iface_num,
    uint8_t f_max,
    uint32_t fwd_table_size,
    __float080 * iface_entry_proportion,
    __float080 * f_distribution) {

    // initialize the network parameters
    this->access_tree_num = access_tree_num;
    this->access_tree_height = access_tree_height;
    this->iface_num = iface_num;
    this->f_max = f_max;
    this->fwd_table_size = fwd_table_size;
    this->iface_entry_proportion = iface_entry_proportion;
    this->f_distribution = f_distribution;

    // initialize the forward() parameters
    this->request_size = 0;
    this->tp_sizes = NULL;
    this->f_r_distribution = NULL;

    // initialize the array of tree roots
    this->root_routers = (RID_Router **) calloc(this->access_tree_num, sizeof(RID_Router *));

    // ... and build the network
    this->build_network();    
}

int RID_Analytics::erase_access_tree_rec(RID_Router * router) {

    uint8_t _iface_num = router->get_iface_num();
    uint8_t _height = router->get_height();

    // recursion ends if we get to the bottom of the access tree
    if (_height >= (this->access_tree_height - 1)) {

        delete router;

        return 0;
    }

    for (int iface = IFACE_DOWNSTREAM; iface < _iface_num; iface++) {

        erase_access_tree_rec(router->get_fwd_table_next_hop(iface));
    }

    delete router;

    return 0;
}

RID_Analytics::~RID_Analytics() {

    // basically erase the path state tree
    tree<Path_State *>::post_order_iterator post_itr = this->path_state_tree.begin_post();
    tree<Path_State *>::post_order_iterator end_post_itr = this->path_state_tree.end_post();

    while (post_itr != end_post_itr) {

        // delete the Path_State object associated with the tree node (not the 
        // node itself)
        delete (*post_itr);

        ++post_itr;
    }

    for (int i = 0; i < this->access_tree_num; i++)
        erase_access_tree_rec(this->root_routers[i]);  

    free(this->root_routers); 
}

int RID_Analytics::print_tp_sizes() {

    if (this->tp_sizes == NULL) {

        fprintf(
            stderr, 
            "RID_Analytics::print_tp_sizes() : tp_sizes is NULL.\n");

        return -1;
    }

    uint32_t _offset_h = 0, _offset_w = 0;

    for (int h = 0; h < this->access_tree_height; h++) {

        _offset_h = h * (pow(2.0, this->access_tree_height - 1) * this->iface_num);

        printf("\nTP SIZES P/ LAYER %d\n", h);

        // top line
        printf("\n-------------");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf("-------------");

        printf("\n[ I ] :     |");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf(" %-11d|", i);

        printf("\n-------------");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf("-------------");

        for (uint8_t w = 0; w < ((int) pow(2.0, h)); w++) {

            printf("\n WIDTH = %-3d|", w);

            _offset_w = w * this->iface_num;

            for (uint8_t i = 0; i < (this->iface_num); i++) {

                printf(" %-11d|", this->tp_sizes[_offset_h + _offset_w + i]);
            }
        }

        printf("\n-------------");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf("-------------");

        printf("\n");
    }

    return 0;
}

int RID_Analytics::print_iface_entry_proportions() {

    if (this->iface_entry_proportion == NULL) {

        fprintf(
            stderr, 
            "RID_Analytics::print_iface_entry_proportions() : iface_entry_proportion is NULL.\n");

        return -1;
    }

    uint32_t _offset_h = 0, _offset_w = 0;

    for (int h = 0; h < this->access_tree_height; h++) {

        _offset_h = h * (pow(2.0, this->access_tree_height - 1) * this->iface_num);

        printf("\nFWD ENTRY PROPORTIONS PER IFACE, LAYER %d\n", h);

        // top line
        printf("\n-------------");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf("-------------");

        printf("\n[ I ] :     |");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf(" %-11d|", i);

        printf("\n-------------");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf("-------------");

        for (uint8_t w = 0; w < ((int) pow(2.0, h)); w++) {

            printf("\n WIDTH = %-3d|", w);

            _offset_w = w * this->iface_num;

            for (uint8_t i = 0; i < (this->iface_num); i++) {

                printf(" %-.5LE|", this->iface_entry_proportion[_offset_h + _offset_w + i]);
            }
        }

        printf("\n-------------");
        for (uint8_t i = 0; i < (this->iface_num); i++)
            printf("-------------");

        printf("\n");
    }

    return 0;
}

int RID_Analytics::build_network_rec(RID_Router * parent_router) {

    uint8_t access_tree = parent_router->get_access_tree_index();
    uint8_t prev_height = parent_router->get_height();
    uint8_t prev_width = parent_router->get_width();

    // recursion ends if we get to the bottom of the access tree
    if (prev_height >= (this->access_tree_height - 1)) {

        // set the bottom-left router of access tree 0 as the starting point 
        // for an RID analytics run
        if (access_tree == 0 && prev_width == 0)
            this->start_router = parent_router;

        return 0;
    }

    uint32_t _offset_h = 0;
    uint32_t _offset_w = 0;

    uint8_t _iface_num = this->iface_num;
    // 'leaf' routers only have IFACE_LOCAL and IFACE_UPSTREAM ifaces. others 
    // will have a full set of ifaces.
    if ((prev_height + 1) == (this->access_tree_height - 1))
        _iface_num = 2;

    // create (this->iface_num - 2) routers below, i.e. for each iface 
    // except IFACE_LOCAL and IFACE_UPSTREAM
    for (uint8_t o = 0; o < (this->iface_num - 2); o++) {

        RID_Router * child_router = 
            new RID_Router(
                access_tree, 
                prev_height + 1,
                (prev_width * (this->iface_num - 2)) + o,
                this->fwd_table_size,
                _iface_num, 
                this->f_max);

        _offset_h = child_router->get_height() * (pow(2.0, this->access_tree_height - 1) * this->iface_num);
        _offset_w = child_router->get_width() * this->iface_num;

        // initialize the ifaces
        for (uint8_t i = 0; i < _iface_num; i++) {
            
            child_router->add_fwd_table_entry(
                i, 
                this->iface_entry_proportion[_offset_h + _offset_w + i], 
                this->f_distribution);

            // printf("RID_Analytics::build_network_rec() : added fwd entry to r[%d][%d] (offset [%d] ?):"\
            //     "\n\t[iface] = %d"\
            //     "\n\t[iface_proportion] = %-.5LE\n", 
            //     child_router->get_height(),  child_router->get_width(), 
            //     _offset_h + _offset_w + i,
            //     i, child_router->get_fwd_table()[i].iface_proportion);
        }

        // we can now set the DOWNSTREAM next hops of the parent router
        parent_router->set_fwd_table_next_hop(IFACE_DOWNSTREAM + o, child_router);
        // set the upstream next hop of the child router
        child_router->set_fwd_table_next_hop(IFACE_UPSTREAM, parent_router);

        build_network_rec(child_router);
    }

    return 0;
}

int RID_Analytics::build_network() {

    for (uint8_t t = 0; t < this->access_tree_num; t++) {

        RID_Router * root_router = 
            new RID_Router(
                t, 0, 0, 
                this->fwd_table_size,
                this->iface_num, 
                this->f_max);

        // initialize the ifaces
        for (uint8_t i = 0; i < this->iface_num; i++) {
            
            root_router->add_fwd_table_entry(
                i, 
                this->iface_entry_proportion[i], 
                this->f_distribution);

            // printf("RID_Analytics::build_network() : added fwd entry to r[%d][%d] (offset [%d] ?):"\
            //     "\n\t[iface] = %d"\
            //     "\n\t[iface_proportion] = %-.5LE\n", 
            //     root_router->get_height(),  root_router->get_width(), 
            //     i,
            //     i, root_router->get_fwd_table()[i].iface_proportion);
        }

        // for now, set the IFACE_UPSTREAM of a root router to NULL 
        root_router->set_fwd_table_next_hop(IFACE_UPSTREAM, NULL);

        // add it to the list of root routers
        this->root_routers[t] = root_router;

        build_network_rec(root_router);
    }

    return 0;
}

int RID_Analytics::run_rec(
    RID_Router * router,
    uint8_t ingress_iface,
    tree<Path_State *>::iterator prev_path_state_itr) {

    // the +2 counts with the special IFACE_MULTI_HITS and IFACE_NO_HITS
    uint8_t _iface_num = router->get_iface_num();
    uint8_t _height = router->get_height();
    uint8_t _width = router->get_width();
    uint8_t _ingress_iface = 0;

    uint32_t _tp_size_offset = 
        _height * (((int) pow(2.0, this->access_tree_height - 1)) * this->iface_num) 
        + _width * this->iface_num;

    __float080 iface_prob = 0.0;
    int iface_path_length = 0;

    // next hop router pointer
    RID_Router * router_next_hop = NULL;
    tree<Path_State *>::iterator path_state_itr;

    // run the forward() operation on this router : this will determine the 
    // possible outgoing ifaces
    printf("\nRID_Analytics::run_rec() : FORWARD() BY ROUTER[%d][%d]\n", _height, _width);

    router->forward(
        this->request_size, 
        ingress_iface,
        &(this->tp_sizes[_tp_size_offset]), 
        this->f_r_distribution);

    // cycle through each of router's ifaces to fill path_state accordingly 
    // and invoke the next hop's forward() (if any router is attached)
    for (int iface = IFACE_MULTI_HITS; iface < _iface_num; iface++) {

        iface_prob = (*prev_path_state_itr)->get_ingress_prob() * router->get_lpm_iface_pmf(iface);
        iface_path_length = (*prev_path_state_itr)->get_path_length() + 1;

        // don't add state nor forward if there's no chance of following this 
        // iface.
        if (iface_prob == 0.0) {
            continue;
        }

        // new path state to be added to the path state tree
        Path_State * path_state = new Path_State(router, iface_prob, iface_path_length);
        path_state_itr = this->path_state_tree.append_child(prev_path_state_itr, path_state);

        // all forwarding outcomes (except IFACE_UPSTREAM and IFACE_UPSTREAM) 
        // result in an end of path (with diff. path outcomes)
        if (iface == IFACE_MULTI_HITS) {

            path_state->set_eop();
            path_state->set_outcome(OUTCOME_MULTI_HITS); 

        } else if (iface == IFACE_NO_HITS) {

            path_state->set_eop();
            path_state->set_outcome(OUTCOME_NO_HITS); 

        } else {

            // a match is verified : it's either have a TP or FP. this 
            // will affect the final outcome if iface == IFACE_LOCAL
            if (this->tp_sizes[_tp_size_offset + iface] > 0) {

                path_state->set_outcome(OUTCOME_TP); 

            } else {

                path_state->set_outcome(OUTCOME_FP);
            }

            if (iface == IFACE_LOCAL) {

                path_state->set_eop();

            } else {

                // get the next hop router by iface
                router_next_hop = router->get_fwd_table_next_hop(iface);

                // this can still happen (or can it?)
                if (router_next_hop == NULL)
                    continue;

                // recursively call run_rec() for router_next_hop : its ingress 
                // iface depends on 'iface' (i.e. the egress iface of 'router'):
                //  * if iface == IFACE_UPSTREAM : ingress iface will be 
                //      IFACE_DOWNSTREAM
                //  * else (i.e. if iface > IFACE_UPSTREAM) : ingress iface 
                //      will be IFACE_UPSTREAM
                _ingress_iface = ((iface == IFACE_UPSTREAM) ? IFACE_DOWNSTREAM : IFACE_UPSTREAM);

                run_rec(router_next_hop, _ingress_iface, path_state_itr);
            }
        } 
    }

    return 0;
}

int RID_Analytics::run(
    uint8_t request_size,
    int * tp_sizes,
    __float080 * f_r_distribution) {

    this->request_size = request_size;
    this->tp_sizes = tp_sizes;
    this->f_r_distribution = f_r_distribution;

    tree<Path_State *>::iterator path_state_tree_itr;

    path_state_tree_itr = this->path_state_tree.begin();
    // FIXME: not sure if it is ok to set the rid_router arg on Path_State() 
    // as NULL
    path_state_tree_itr = this->path_state_tree.insert(
                                                    path_state_tree_itr, 
                                                    new Path_State(NULL, 1.0, 0));

    // ingress iface is the IFACE_DOWNSTREAM of upstream router
    run_rec(this->start_router, IFACE_DOWNSTREAM, path_state_tree_itr);

    return 0;
}

int RID_Analytics::view_results(
    uint8_t modes,
    char * output_file_path) {

    // iterate through the path_tree leafs to learn 
    // a bit about latencies for this scenario
    tree<Path_State *>::leaf_iterator leaf_itr = this->path_state_tree.begin_leaf();
    tree<Path_State *>::leaf_iterator end_leaf_itr = this->path_state_tree.end_leaf();

    __float080 checksum = 0.0;

    char * path_state_str = NULL;

    // print some results for quick checking
    if ((modes & MODE_VERBOSE))
        printf("\n*** RESULTS ***\n\n");

    // if requested, save outcomes in the path given as arg
    FILE * output_file;

    if ((modes & MODE_SAVEOUTCOMES)) {

        output_file = fopen(output_file_path, "a");
    }

    while (leaf_itr != end_leaf_itr) {

        if ((*leaf_itr)->get_ingress_prob() > 0.0) {

            if ((modes & MODE_VERBOSE)) {

                path_state_str = (*leaf_itr)->to_string();

                printf("[PATH_TREE DEPTH %d] : %s\n", 
                    this->path_state_tree.depth(leaf_itr), 
                    path_state_str);

                free(path_state_str);
            }

            if ((modes & MODE_SAVEOUTCOMES)) {

                // append a line in the output file, with the following 
                // format : 
                // <rtr.hght>,<rtr.wdth>,<prob>,<path_length>,<outcome>
                fprintf(output_file, 
                    "%d,%d,%-.5LE,%d,%d\n", 
                    (*leaf_itr)->get_router()->get_height(),
                    (*leaf_itr)->get_router()->get_width(),
                    (*leaf_itr)->get_ingress_prob(),
                    (*leaf_itr)->get_path_length(),
                    (*leaf_itr)->get_outcome());
            }

            checksum += (*leaf_itr)->get_ingress_prob();
        }

        ++leaf_itr;
    }

    if ((modes & MODE_VERBOSE)) {

        printf("\n[CHECKSUM : %-.5LE]\n", checksum);
    }

    return 0;
}
