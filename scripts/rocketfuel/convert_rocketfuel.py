import argparse
import re
import sys
import glob
import os
import math
import binascii
import random

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import xml.etree.cElementTree as et

from xml.dom import minidom
from prettytable import PrettyTable

from collections import defaultdict
from collections import OrderedDict
from random import randint

DEFAULT_DIFF_DISTR_DIR = "/home/adamiaonr/workbench/rid-analytics/experiments/examples/diff-distributions"

# FIXME : this algorithm is terribly inefficient...
def get_fwd_dist(topology, shortest_paths, router):

    # fwd entry size distr. values
    fwd_dist = defaultdict(float)
    # holds the nodes accessible via shortest paths, through each neighbor 
    iface_trees = defaultdict(set)

    # nr. of nodes in topology
    node_num = len(topology.nodes())
    neighbor_num = len(topology.neighbors(router))
    # set of visited sources
    visited_sources = set([])

    # how many nodes announce their contents via neighbor?
    # we determine this by finding the nr. of sources which include a shortest 
    # paths with the sequence "neighbor,router" in it.
    for i, neighbor in enumerate(topology.neighbors(router)):
        
        search_str = ("%03d,%03d" % (neighbor, router))
        
        for source in shortest_paths:

            if source in visited_sources:
                continue

            for path in shortest_paths[source]:

                if search_str in path:
                    fwd_dist[neighbor] += 1.00
                    visited_sources.add(source)
                    iface_trees[neighbor].add(source)
                    break

    for neighbor in fwd_dist:
        fwd_dist[neighbor] = fwd_dist[neighbor] * (1.0 / float(node_num))

    # interface 0 (local) has 1.0 / neighbor_num fwd_dist value
    fwd_dist[router] = (1.0 / float(node_num))

    return fwd_dist, iface_trees

def to_tree_bitmask(topology, iface_trees, entry_sizes):

    # calculate nr. of bytes of the bitmask. 
    # 8 nodes per byte
    nr_bytes = int(math.ceil(float(len(topology.topology.nodes())) / 8.0))

    # initialize tree bitmasks
    # we use 2, for convience (i.e. it's more intuitive to do 
    # this in python than c++):
    #   - tree_bitmasks[s], per each fwd. entry size
    #   - tree_bitmask, union of all sizes
    tree_bitmasks = defaultdict(list)
    tree_bitmasks[0] = [0] * nr_bytes

    # 1) start w/ base sizes (i.e. info in entry sizes, which 
    # is applicable to all nodes)
    for s in entry_sizes:

        if not (entry_sizes[s] > 0.0):
            continue

        tree_bitmasks[s] = [0] * nr_bytes
        for source in iface_trees:
            byte_index = int((source / 8))
            # print("source %d, byte index : %d" % (source, byte_index))
            # print("source : %d, tree_bitmask[%d] = %d" % (source, byte_index, tree_bitmask[byte_index]))
            tree_bitmasks[s][byte_index] |= (1 << (source % 8))
            # print("bit index = %d, tree_bitmask[%d] = %d" % ((source % 8), byte_index, tree_bitmask[byte_index]))
            tree_bitmasks[0][byte_index] |= (1 << (source % 8))

    # 2) add fp route contributions
    for s in topology.fp_srcs:

        srcs = topology.fp_srcs[s].intersection(set([src for src in iface_trees]))

        if s not in tree_bitmasks:
            tree_bitmasks[s] = [0] * nr_bytes

        for src in srcs:
            byte_index = int((src / 8))
            tree_bitmasks[s][byte_index] |= (1 << (src % 8))
            tree_bitmasks[0][byte_index] |= (1 << (src % 8))

    # 3) add tp route contributions

    # 4) 'hexlify' the treebitmasks
    for s in tree_bitmasks:
        tree_bitmasks[s] = binascii.hexlify(bytearray(tree_bitmasks[s]))

    return nr_bytes, tree_bitmasks

def load_request_entry_diff_distributions(req_size, topology_block, diff_distr_type = 'l', diff_distr_dir = DEFAULT_DIFF_DISTR_DIR):

    filename = os.path.join(diff_distr_dir, ("%d.fdist" % (req_size)))
    for line in open(filename, "r").readlines():

        if line.split("=", 1)[0] is not diff_distr_type:
            continue

        values = line.split("=", 1)[1].split(",")
        for i, value in enumerate(values):
            if (i + 1) > req_size:
                break
            f_r_dist_block = et.SubElement(topology_block, "f_r_dist", diff = str(i + 1)).text = value

# post simple statistics about each AS-level topology:
#   - nr. of nodes
#   - nr. of links
#   - median outdegree
def print_pop_level_statistics(rocketfuel_dir):

    stats = OrderedDict()

    for root, dirs, files in os.walk(rocketfuel_dir):
        for filename in files:

            if not filename.endswith(".wt"):
                continue

            filename = os.path.join(root, filename)
            print("filename : %s" % (filename))

            topology_id = filename.split("/")[-2]
            topology_id = topology_id.split(".", 1)[0]

            topology_obj = parse_pop_level_map(filename)
            stats[topology_id] = topology_obj.get_pop_level_statistics()

    table = PrettyTable(['as id', '# nodes', '# links', 'min. outdegree', 'median outdegree', 'mean outdegree', 'max. outdegree', 'median length', 'max. length'])

    for topology_id in stats:
        table.add_row([
            topology_id,
            stats[topology_id]['nodes'],
            stats[topology_id]['links'],
            stats[topology_id]['min-outdegree'],
            stats[topology_id]['med-outdegree'],
            stats[topology_id]['mean-outdegree'],
            stats[topology_id]['max-outdegree'],
            (np.median(stats[topology_id]['path-lengths']) - 1),
            (stats[topology_id]['path-lengths'][0] - 1)
            ])

    print("")
    print(table)
    print("")

class Topology:

    def __init__(self, topology):

        self.topology = topology
        self.shortest_paths = defaultdict(list)
        self.stats = defaultdict()

        # tp and fp sources - added w/ add_fp_route() and add_tp_route() - 
        # indexed by size
        self.tp_srcs = defaultdict(set)
        self.fp_srcs = defaultdict(set)

    def draw_pop_level_map(self, filename, filters = {}):

        # setting the positions with respect to the whole graph
        # (not just a subgraph)
        pos = nx.spring_layout(self.topology)

        labels = {}
        for router in self.topology.nodes():
            labels[router] = str(router)

        if 'annc-size' in filters:
            # collect nodes by annc. size they send / recv
            annc_sizes = {'fp' : defaultdict(set), 'tp' : defaultdict(set)}
            for router in self.topology.nodes():
                for i, neighbor_router in enumerate(self.topology.neighbors(router)):

                    # print("%s <->> %s" % (router, neighbor_router))
                    
                    edge = self.topology.get_edge_data(router, neighbor_router)
                    if 'fp' in edge:
                        for fp_record in edge['fp'].split("|"):
                            # print("fp : %s" % (fp_record))
                            annc_sizes['fp'][int(fp_record.split(':')[2])].add(router)
                            annc_sizes['fp'][int(fp_record.split(':')[2])].add(neighbor_router)

                    if 'tp' in edge: 
                        for tp_record in edge['tp'].split("|"):
                            # print("tp : %s" % (tp_record))
                            annc_sizes['tp'][int(tp_record.split(':')[2])].add(router)
                            annc_sizes['tp'][int(tp_record.split(':')[2])].add(neighbor_router)

            colors = {1 : 'black', 2 : 'red', 4 : 'red', 5 : 'blue', 10 : 'green'}
            for s in annc_sizes['fp']:
                
                # generate subgraph
                sg = self.topology.subgraph(list(annc_sizes['fp'][s]))

                nx.draw_networkx_nodes(
                    sg, 
                    pos = pos,
                    node_color = 'red',
                    node_size = 500,
                    alpha = 0.5)

                nx.draw_networkx_edges(
                    sg,
                    pos,
                    edge_color = colors[s],
                    width = 1.0, 
                    alpha = 0.5)

                nx.draw_networkx_labels(sg, pos = pos, labels = labels, font_size = 10)

            # for s in annc_sizes['tp']:
                
            #     # generate subgraph
            #     sg = self.topology.subgraph(list(annc_sizes['tp'][s]))

            #     nx.draw_networkx_nodes(
            #         sg, 
            #         pos = pos,
            #         node_color = 'red',
            #         node_size = 500,
            #         alpha = 0.5)

            #     nx.draw_networkx_edges(
            #         sg,
            #         pos,
            #         edge_color = colors[s],
            #         edge_style = 'dashed',
            #         width = 1.0, 
            #         alpha = 0.5)

            #     nx.draw_networkx_labels(sg, pos = pos, labels = labels, font_size = 10)

        # nx.draw_networkx_nodes(
        #     self.topology, 
        #     pos,
        #     node_color = 'red',
        #     node_size = 500,
        #     alpha = 0.5)

        # nx.draw_networkx_edges(
        #     self.topology,
        #     pos,
        #     width = 1.0, 
        #     alpha = 0.5)

        # labels = {}
        # for router in self.topology.nodes():
        #     labels[router] = str(router)
        # nx.draw_networkx_labels(self.topology, pos, labels, font_size = 10)

        # save the figure in <rocketfuel-file>.pdf
        plt.savefig(filename, bbox_inches='tight', format = 'pdf')

    def get_shortest_paths(self, force = False):

        if (len(self.shortest_paths) > 0) and (force == False):
            return self.shortest_paths

        dsts = defaultdict(set)

        # compute shortest paths from each node to every other node.
        # paths are stored as strings, sequences of node ids, separated by ','. 
        # the structure of the returned dictionary is:
        #
        #   - [src 1] -> [ 
        #                   "<src 1>,<path node 1>,<path node 2>, ... ,<node 0>",
        #                   "<src 1>,<path node 1>,<path node 2>, ... ,<node 1>",
        #                   (...),
        #                   "<src 1>,<path node 1>, ... ,<node n>"
        #                   ]
        #
        #   - [src 2] -> [ 
        #                   "<src 2>,<path node 1>,<path node 2>, ... ,<node 0>",
        #                   "<src 2>,<path node 1>,<path node 2>, ... ,<node 1>",
        #                   (...),
        #                   "<src 2>,<path node 1>, ... ,<node n>"
        #                   ]
        #   - (...)
        #
        self.shortest_paths = defaultdict(list)

        for source in self.topology.nodes():

            routes = nx.single_source_shortest_path(self.topology, source)

            if source not in dsts:
                dsts[source] = set([])

            for dst, route in routes.iteritems():

                if not ((source in self.shortest_paths) and (dst in dsts[source])):
                    self.shortest_paths[source].append(','.join(("%03d" % (node)) for node in route))
                    dsts[source].add(dst)

                # we assume routes are symmetric :
                #   - i.e., the shortest path between a -> b, is the same as b -> a
                #   - to force it, we add the route from dst to source as 
                #     the inverted _route
                if dst not in dsts:
                    dsts[dst] = set([])

                if not ((dst in self.shortest_paths) and (source in dsts[dst])):
                    self.shortest_paths[dst].append(','.join(("%03d" % (node)) for node in reversed(route)))
                    dsts[dst].add(source)

        return self.shortest_paths

    def set_shortest_paths(self, shortest_paths):
        self.shortest_paths = shortest_paths

    def get_neighbors_within_dist(self, node_id, dist):
        return [node for node in nx.single_source_shortest_path(self.topology, node_id, dist)]

    # add true positive entries to a topology (.scn file). 
    # we start by applying changes in a networkx object, then call 
    # convert_to_scn() to generate a .scn file.
    #
    # algorithm:
    #   - find all shortest paths in topology originated on tp_src_id 
    #   - add entries of size tp_size on links in shortest path
    def add_tp_route(self, tp_src_id, tp_size = 1, radius = -1, n_srcs = 1, annc_to = -1, avoid_path = []):

        # if tp_src_id isn't defined, we pick n_srcs among the radius-hop neighbors of annc_to
        if tp_src_id == -1:

            # extract neighbors of annc_to at radius hops from it
            neighbors = self.get_neighbors_within_dist(annc_to, radius)
            # select n_srcs neighbors
            for neighbor in neighbors:

                if n_srcs < 1:
                    return 0

                # find the shortest paths from neighbor to annc_to by inspecting
                # the shortest paths which:
                #   - depart from neighbor
                #   - end in annc_to
                #
                # this is guaranteed to be the shortest path from annc_to to neighbor, since
                # shortest paths are symmetric, i.e. path(neighbor > annc_to) = path(annc_to > neighbor)
                for path in self.get_shortest_paths()[neighbor]:

                    path = [ int(p) for p in path.split(",") ]

                    # if the end of that path isn't annc_to, continue
                    if path[-1] != annc_to:
                        continue

                    # if neighbor is below radius hops, abort
                    if ((len(path) - 1) < radius):
                        print("%d -> %d : %s (%d : no radius)" % (path[0], path[-1], str(path), len(path) - 1))
                        continue

                    # if neighbor is in avoid_path, continue
                    if neighbor in avoid_path:
                        print("potential tp src in avoid path : %d vs. [%s]. aborting." % (neighbor, avoid_path))
                        continue

                    print("adding secondary tp source : %d -> %d : %s (size : %d)" % (path[0], path[-1], str(path), tp_size))

                    # add path[0] - i.e. the origin of the fp route - to the fp_srcs dict
                    self.tp_srcs[tp_size].add(path[0])

                    # add a tp entry to the local iface of tp_src_id. this is basically 
                    # an edge with the same head and tail.

                    # preserve the tp attributes on edges
                    edge = self.topology.get_edge_data(neighbor, neighbor)

                    # we will save the pre-existing tp info on a dict
                    edge_tps = defaultdict()
                    if (edge is not None) and ('tp' in edge):
                        for tp_record in edge['tp'].split("|"):
                            edge_tps[':'.join(tp_record.split(":")[:2])] = tp_record.split(":")[2]

                    # prepare a key for the local iface on neighbor
                    edge_tps_key = ("%d:%d" % (neighbor, 0))
                    if edge_tps_key in edge_tps:
                        # replace the existing tp_size value only if x < tp_size
                        x = int(edge_tps[edge_tps_key])
                        if tp_size > x:
                            edge_tps[edge_tps_key] = str(tp_size)
                    else:
                        edge_tps[edge_tps_key] = str(tp_size)

                    # convert edge_tps back to string form
                    edge_tps_str = ""
                    for edge_tps_key in edge_tps:
                        edge_tps_str += ("%s:%s|" % (edge_tps_key, edge_tps[edge_tps_key]))

                    self.topology.add_edge(neighbor, neighbor, 
                        type = 'internal', 
                        e1 = ("%d:%d"       % (neighbor, 0)), 
                        e2 = ("%d:%d"       % (neighbor, 0)),
                        tp = (edge_tps_str.rstrip("|")))

                    # iterate through the links in the path, so that we add tp information to them
                    for i in xrange(len(path) - 1):
                        
                        # get src and dst ends of the edge
                        src, dst = path[i], path[i + 1]
                        # determine the interface of dst to which we should add the tp info
                        edge = self.topology.get_edge_data(src, dst)

                        # we will save the pre-existing tp info on a dict
                        edge_tps = defaultdict()
                        if (edge is not None) and ('tp' in edge):
                            for tp_record in edge['tp'].split("|"):
                                edge_tps[':'.join(tp_record.split(":")[:2])] = tp_record.split(":")[2]

                        iface = 0
                        if int(edge['e1'].split(':', 1)[0]) == dst:
                            iface = edge['e1']
                        else:
                            iface = edge['e2']

                        # prepare a key for the local iface on neighbor
                        edge_tps_key = ("%s" % (iface))
                        # print(edge_tps_key)
                        if edge_tps_key in edge_tps:
                            # replace the existing tp_size value only if x < tp_size
                            x = int(edge_tps[edge_tps_key])
                            if tp_size > x:
                                edge_tps[edge_tps_key] = str(tp_size)
                        else:
                            edge_tps[edge_tps_key] = str(tp_size)

                        # for edge_tps_key in edge_tps:
                        #     print("%s : %s" % (edge_tps_key, edge_tps[edge_tps_key]))

                        # convert edge_tps back to string form
                        edge_tps_str = ""
                        for edge_tps_key in edge_tps:
                            edge_tps_str += ("%s:%s|" % (edge_tps_key, edge_tps[edge_tps_key]))

                        # print("%d -> %d : %s" % (src, dst, edge_tps_str))

                        # the tp info is a tuple of the form <router>:<outgoing iface>:<size> 
                        self.topology.add_edge(src, dst, tp = edge_tps_str.rstrip("|"))

                    n_srcs -= 1
                    break

            return 0

        print("Topology::add_tp_route() : [INFO] adding tp routes towards %s" % (str(tp_src_id)))
        
        # get all shortest path starting at tp_src
        routes = self.get_shortest_paths()[tp_src_id]
        # add a tp entry to the *local* iface of tp_src
        self.topology.add_edge(tp_src_id, tp_src_id, 
            type = 'internal', 
            e1 = ("%d:%d"       % (tp_src_id, 0)), 
            e2 = ("%d:%d"       % (tp_src_id, 0)),
            tp = ("%d:%d:%d"    % (tp_src_id, 0, tp_size)))

        # iterate through the links in the path, add tp information to them
        for route in routes:

            path = [ int(n) for n in route.split(",") ]
            # print("Topology::add_tp_route() : [INFO] adding tp route %s" % (str(path)))

            for i in xrange(len(path) - 1):
                # get src and dst ends of the edge
                src, dst = path[i], path[i + 1]
                 # determine the interface of dst to which we should add the tp info
                edge = self.topology.get_edge_data(src, dst)

                iface = 0
                if int(edge['e1'].split(':', 1)[0]) == dst:
                    iface = edge['e1']
                else:
                    iface = edge['e2']

                # the tp info is a tuple of the form <router>:<outgoing iface>:<size> 
                self.topology.add_edge(src, dst, tp = ("%s:%d" % (iface, tp_size)))

        return 0

    def add_fp_route(self, src_id, fp_size = 2, fp_size_proportion = 10, annc_radius = 1, dist = 1, avoid_path = [], n_srcs = [0,0]):

        # get neighbors within dist of src_id
        neighbors = self.get_neighbors_within_dist(src_id, dist)

        n_neighbors = len(neighbors)
        if (n_srcs[0] == 2):
            n_srcs[1] = n_neighbors

        # iterate over neighbors of src_id
        for neighbor in neighbors:

            # if specified, add a limited number of sources
            if (n_srcs[0] > 0) and (n_srcs[1] < 1):
                return 0

            if (n_srcs[0] == 2) and (n_srcs[1] < (n_neighbors / 2)):
                return 0

            # skip possible fp caches w/ 50% probability if n_srcs[0] is of mode 2
            if (n_srcs[0] > 0) and (bool(random.getrandbits(1))):
                continue

            # iterate over all paths starting in the neighbor
            for path in self.get_shortest_paths()[neighbor]:

                path = [ int(p) for p in path.split(",") ]

                # if neighbor is below dist hops from source, abort
                if ((len(path) - 1) < dist):
                    # print("%d -> %d : %s (%d : no dist)" % (path[0], path[-1], str(path), len(path) - 1))
                    continue

                # # if the end of that path isn't src_id, continue
                # if path[-1] != src_id:
                #     continue

                # if neighbor is in avoid_path, continue
                if neighbor in (avoid_path[:dist] + avoid_path[dist + 1:]):
                    # print("Topology::add_fp_route() : [INFO] potential fp src in avoid path : %d vs. [%s]. aborting." % (neighbor, avoid_path))
                    continue

                # add path[0] - i.e. the origin of the fp route - to the fp_srcs dict
                self.fp_srcs[fp_size].add(path[0])
                # print("Topology::add_fp_route() : [INFO] adding fp route : %d -> %d : %s (%d)" % (path[0], path[-1], str(path), len(path) - 1))

                # update the nr of srcs added so far
                if (n_srcs[0] > 0):
                    n_srcs[1] = (n_srcs[1] - 1)

                # add a fp entry to the local iface of neighbor. 
                # this is basically an edge with the same head and tail
                # get any fp edge attribute, if already existent. 
                # the point is to update it if a previous entry for the same iface:size already exists
                edge = self.topology.get_edge_data(neighbor, neighbor)
                # we will save the pre-existing fp info on a dict
                edge_fps = defaultdict()
                if (edge is not None) and ('fp' in edge):
                    for fp_record in edge['fp'].split("|"):
                        edge_fps[':'.join(fp_record.split(":")[:3])] = fp_record.split(":")[3]

                # prepare a key for the local iface on neighbor
                edge_fps_key = ("%d:%d:%d" % (neighbor, 0, fp_size))
                if edge_fps_key in edge_fps:
                    x = int(edge_fps[edge_fps_key])
                    edge_fps[edge_fps_key] = str(x + fp_size_proportion)
                else:
                    edge_fps[edge_fps_key] = str(fp_size_proportion)

                # convert edge_fps back to string form
                edge_fps_str = ""
                for edge_fps_key in edge_fps:
                    edge_fps_str += ("%s:%s|" % (edge_fps_key, edge_fps[edge_fps_key]))

                # add an update 'fp' attribute to the local iface edge 
                # print(edge_fps_str)
                self.topology.add_edge(neighbor, neighbor, 
                    e1 = ("%d:%d"   % (neighbor, 0)), 
                    e2 = ("%d:%d"   % (neighbor, 0)),
                    fp = (edge_fps_str.rstrip("|")))

                # print("<router id='%d'>" % (int(neighbor)))
                # print("add <fwd_size_dist size='%d'>%.2f</fwd_size_dist> to <link local='%d'>\n" 
                #     % (fp_size, (float(fp_size_proportion) / 100.0), 0))

                # for i in xrange(len(path) - 1):
                for i in xrange(annc_radius):

                    # make sure we don't over the end of a path
                    if (i + 1) > (len(path) - 1):
                        break

                    # get src and dst ends of the edge. 
                    # the announcement works from src -> dst, as such, it should be stored in the iface of dst.
                    src, dst = path[i], path[i + 1]

                    # determine the interface of dst to which we should add the fp info
                    edge = self.topology.get_edge_data(src, dst)
                    # retrieve any pre-existing fp information in the form of a dict(), if any
                    edge_fps = defaultdict()
                    if 'fp' in edge:
                        # print("edge[fp] = %s" % (edge['fp']))
                        for fp_record in edge['fp'].split("|"):
                            # print("fp_record = %s" % (fp_record))
                            # print("manufactured key = %s" % (':'.join(fp_record.split(":")[:3])))
                            edge_fps[':'.join(fp_record.split(":")[:3])] = fp_record.split(":")[3]

                    # find out the iface nr. on dst node, to which the fp info should be added to
                    iface = 0
                    if int(edge['e1'].split(':', 1)[0]) == dst:
                        iface = edge['e1']
                    else:
                        iface = edge['e2']

                    # if a fp record for fp_size already exists on this iface, 
                    # update the proportion value. if not, just add a new key
                    edge_fps_key = ("%s:%d" % (iface, fp_size))
                    
                    # print("key = %s" % (edge_fps_key))
                    # print("keys = %s" % (str([k for k in edge_fps])))

                    if edge_fps_key in edge_fps:
                        x = int(edge_fps[edge_fps_key])
                        edge_fps[edge_fps_key] = str(x + fp_size_proportion)
                    else:
                        edge_fps[edge_fps_key] = str(fp_size_proportion)
                    # convert edge_fps back to string form
                    edge_fps_str = ""
                    for edge_fps_key in edge_fps:
                        edge_fps_str += ("%s:%s|" % (edge_fps_key, edge_fps[edge_fps_key]))
                    # the tp info is a tuple of the form <router>:<outgoing iface>:<size> 

                    self.topology.add_edge(src, dst, fp = edge_fps_str.rstrip("|"))
                    # print("<router id='%d'>" % (int(dst)))
                    # print("add <fwd_size_dist size='%d'>%.2f</fwd_size_dist> to <link local='%d'>\n" 
                    #     % (fp_size, (float(fp_size_proportion) / 100.0), int(iface.split(":")[1])))

        # for edge in self.topology.edges():
        #     print(self.topology.get_edge_data(edge[0], edge[1]))

        return 0

    def generate_paths(self, nr_paths, path_size):

        stats = self.get_pop_level_statistics()
        shortest_paths = self.get_shortest_paths()

        print("Topology::generate_paths() : [INFO] %d shortest path sources for [n : %d, size : %d]:" 
            % (len(shortest_paths), nr_paths, path_size))

        # collect x valid paths per source
        paths = defaultdict()
        sources = set([])
        for source in reversed([s for s in shortest_paths]):

            n_paths = 2
            for path in [p.split(",") for p in shortest_paths[source]]:

                curr_path_size = (len(path) - 1)

                if curr_path_size == path_size:
                    sources.add(source)
                    avg_path_outdegree = float(sum([stats['outdegree-list'][int(r)] for r in path])) / float(curr_path_size)
                    paths[("%d:%d:%.2f" % (int(path[0]), int(path[-1]), avg_path_outdegree))] = [int(p) for p in path]

                    n_paths -= 1
                    if (n_paths < 1):
                        break

        print("Topology::generate_paths() : [INFO] collected %d paths from %d diff. sources" 
            % (len(paths), len(sources)))

        # pick nr_paths ar random from paths
        keys = [k for k in paths]
        random.shuffle(keys)
        return defaultdict(list, ((k, paths[k]) for k in keys[0:nr_paths]))

    def get_paths(self, selected_paths):

        paths = OrderedDict()
        stats = self.get_pop_level_statistics()
        shortest_paths = self.get_shortest_paths()

        for selected_path in selected_paths:
            for path in [p.split(",") for p in shortest_paths[int(selected_path.split(":")[0])]]:
                # stop when path[-1] is == selected_path[-1]
                if int(path[-1]) == int(selected_path.split(":")[-1]):
                    avg_path_outdegree = float(sum([stats['outdegree-list'][int(r)] for r in path])) / float(len(path) - 1)
                    paths[("%d:%d:%.2f" % (int(path[0]), int(path[-1]), avg_path_outdegree))] = [int(p) for p in path]
                    break

        print("Topology::get_paths() : [INFO] collected %d paths, as requested:" % (len(paths)))
        for i, path in enumerate(paths):
            print("\t[%d:%d] : %s" % (
                int(selected_paths[i].split(":")[0]), int(selected_paths[i].split(":")[-1]), 
                str(paths[path])))

        return paths

    def get_pop_level_statistics(self, force = False):

        if (len(self.stats) > 0) and (force == False):
            return self.stats

        self.stats = defaultdict()

        self.stats['nodes'] = self.topology.number_of_nodes()
        self.stats['links'] = self.topology.number_of_edges()
        self.stats['outdegree-list'] = nx.degree(self.topology)
        degree_sequence = sorted(nx.degree(self.topology).values(), reverse = True)

        avg_path_outdegrees = []
        path_lengths = []
        _shortest_paths = self.get_shortest_paths()
        for source in _shortest_paths:
            for path in _shortest_paths[source]:

                path_lengths.append(len(path.split(",")))
                path_routers = [int(r) for r in path.split(",")]
                avg_path_outdegrees.append(float(sum([self.stats['outdegree-list'][r] for r in path_routers])) / (float(len(path.split(",")))))

        self.stats['avg-path-outdegree'] = list(set(avg_path_outdegrees))
        self.stats['path-lengths'] = sorted(list(set(path_lengths)), reverse = True)
        self.stats['min-outdegree'] = np.min(degree_sequence)
        self.stats['med-outdegree'] = np.median(degree_sequence)
        self.stats['mean-outdegree'] = np.mean(degree_sequence)
        self.stats['max-outdegree'] = np.max(degree_sequence)

        return self.stats


    def convert_to_scn(self, entry_sizes, table_size = 100000000, req_size = 5, scn_filename = 'topology.scn'):

        # add <topology> block (main block of .scn file)
        topology_block = et.Element("topology")

        # get shortest path list
        shortest_paths = self.get_shortest_paths()
        # determine a ttl value for the topology 
        # 2 * longest path, as in http://www.map.meteoswiss.ch/map-doc/ftp-probleme.htm
        ttl = 0
        for source in shortest_paths:
            for path in shortest_paths[source]:
                path_length = (len(path.split(",")) - 1)
                if ttl < path_length:
                    ttl = path_length

        # add <ttl> block
        ttl_block = et.SubElement(topology_block, "ttl").text = str(ttl)
        # load |F\R| distributions from .fdist files & add <f_r_dist diff> blocks,
        # one per possible |F\R| value
        load_request_entry_diff_distributions(req_size, topology_block)
        # add <fwd_table_size> block
        fwd_table_size_block = et.SubElement(topology_block, "fwd_table_size").text = str(table_size)

        # cycle through each node in the topology. 
        # we then write 5 diff. types of blocks: 
        #   - <router id=""> block
        #   - <link local="" remote="" rrouter=""> block
        #   - <fwd_size_dist size="">x.xx</fwd_size_dist> block
        #   - <fwd_dist iface="y">x.xx</fwd_dist>
        #   - <tree_bitmask len="" size=""> ... </tree_bitmask>
        #   - <tp iface="y">x</tp>
        for router in self.topology.nodes():

            added_local_iface = False

            # add <router id=""> block
            router_block = et.SubElement(topology_block, "router", id = str(router))

            # get: 
            #   - the fwd entry size distr. values for each of the ifaces in router
            #   - the set of nodes accessible via each of the router's neighbors
            fwd_dist, iface_trees = get_fwd_dist(self.topology, shortest_paths, router)

            # add <link> blocks, one per each neighbor of router
            for i, neighbor_router in enumerate(self.topology.neighbors(router)):

                # extract the attributes of the link
                edge = self.topology.get_edge_data(router, neighbor_router)

                # distribute the link ifaces 'a' and 'b' over local and remote
                local_iface     = 0
                remote_iface    = 0

                if int(edge['e1'].split(':', 1)[0]) == router:
                    local_iface   = edge['e1'].split(':', 1)[1]
                    remote_iface  = edge['e2'].split(':', 1)[1]
                else:
                    local_iface   = edge['e2'].split(':', 1)[1]
                    remote_iface  = edge['e1'].split(':', 1)[1]

                link_block = et.SubElement(router_block, "link", 
                    local   = str(local_iface), 
                    remote  = str(remote_iface), 
                    rrouter = str(neighbor_router))

                # add <fwd_size_dist> blocks : distribution of forwarding entry sizes at each link
                fwd_size_dist = defaultdict(float)
                fwd_size_dist_total = 0.0

                # check for 'fp' entries in the link
                if 'fp' in edge:
                    fp_records = edge['fp'].split("|")
                    for fp_record in fp_records:
                        if ((int(fp_record.split(':')[0]) == router) and (int(fp_record.split(':')[1]) == int(local_iface))):
                            fwd_size_dist[int(fp_record.split(':')[2])] += (float(fp_record.split(':')[3]) / 100.0)
                            fwd_size_dist_total += (float(fp_record.split(':')[3]) / 100.0)

                for i in xrange(req_size):
                    if str(i + 1) in entry_sizes:
                        fwd_size_dist[i + 1] += float(entry_sizes[str(i + 1)])
                        fwd_size_dist_total  += float(entry_sizes[str(i + 1)])

                for i in xrange(req_size):
                    fwd_size_dist_block_text = "0.00"
                    if (i + 1) in fwd_size_dist:
                        fwd_size_dist_block_text = ("%.2f" % (fwd_size_dist[i + 1] / fwd_size_dist_total))
                    fwd_size_dist_block = et.SubElement(link_block, "fwd_size_dist", size = str(i + 1)).text = fwd_size_dist_block_text

                # add <tp> block (if edge has 'tp' attribute)
                if ('tp' in edge): 
                    for tp_record in edge['tp'].split("|"):
                        if int(tp_record.split(':')[0]) == router:
                            tp_block = et.SubElement(router_block, "tp", iface = tp_record.split(':')[1]).text = tp_record.split(':')[2]

                # add <tree_bitmask> blocks, one per entry size (only entry sizes for which % > 0)
                if (neighbor_router == router):

                    # adding local iface

                    nr_bytes, tree_bitmasks = to_tree_bitmask(self, [router], entry_sizes)

                    # - tree bitmask per fwd entry size
                    for s in tree_bitmasks:
                        tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", len = str(nr_bytes), size = str(s)).text = tree_bitmasks[s]
    
                    # mark local iface as added
                    added_local_iface = True
                
                else:
                    # adding any other iface
                    nr_bytes, tree_bitmasks = to_tree_bitmask(self, list(iface_trees[neighbor_router]), entry_sizes)

                    for s in tree_bitmasks:
                        tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", len = str(nr_bytes), size = str(s)).text = tree_bitmasks[s]

            # add a special local link (represents local delivery), if not added already
            if not added_local_iface:
                link_block = et.SubElement(router_block, "link", 
                    local   = str(0), 
                    remote  = str(0), 
                    rrouter = str(router))

                for i in xrange(req_size):

                    fwd_size_dist_block_text = "0.00"
                    if str(i + 1) in entry_sizes:
                        fwd_size_dist_block_text = ("%.2f" % (entry_sizes[str(i + 1)]))

                    fwd_size_dist_block = et.SubElement(link_block, "fwd_size_dist", size = str(i + 1)).text = fwd_size_dist_block_text

                nr_bytes, tree_bitmasks = to_tree_bitmask(self, [router], entry_sizes)
                for s in tree_bitmasks:
                    tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", len = str(nr_bytes), size = str(s)).text = tree_bitmasks[s]

            # add <fwd_dist> blocks for each iface
            # if the local iface isn't calculated in get_fwd_dist() (e.g. that 
            # only happens if the link (router, router) has been added), add it now.
            if not added_local_iface:
                fwd_dist_block = et.SubElement(router_block, "fwd_dist", iface = str(0)).text = ("%.4f" % (fwd_dist[router]))

            # finally, the ifaces pointing to each neighboring router
            for neighbor_router in self.topology.neighbors(router):

                # find iface which links to neighbor_router
                edge = self.topology.get_edge_data(router, neighbor_router)
                local_iface = 0
                if int(edge['e1'].split(':', 1)[0]) == router:
                    local_iface = edge['e1'].split(':', 1)[1]
                else:
                    local_iface = edge['e2'].split(':', 1)[1]

                fwd_dist_block = et.SubElement(router_block, "fwd_dist", iface = str(local_iface)).text = ("%.4f" % (fwd_dist[neighbor_router]))


        xmlstr = minidom.parseString(et.tostring(topology_block)).toprettyxml(indent="    ")
        with open(scn_filename, "w") as f:
            f.write(xmlstr)

    def trim_topology(self, max_outdegree):

        # extract stats from current topology
        stats = self.get_pop_level_statistics()

        # check the outdegree list, look for nodes w/ outdegree > max_outdegree
        for router, outdegree in stats['outdegree-list'].iteritems():

            if outdegree <= max_outdegree:
                continue

            # get the nr. of neighbors of router in excess to trim
            to_trim = outdegree - max_outdegree

            # 1st remove neighbors w/ outdegree 1. if there are still neighbors 
            # left, remove those w/ outdegree 2, and so on.
            outdegrees_for_removal = [1, 2, 3]
            for r in outdegrees_for_removal: 

                if to_trim == 0:
                    break

                for i, neighbor_router in enumerate(self.topology.neighbors(router)):

                    if to_trim == 0:
                        break

                    if stats['outdegree-list'][neighbor_router] == r:
                        self.topology.remove_node(neighbor_router)
                        to_trim -= 1

        # remove nodes with no neighbors
        for node in self.topology.nodes():

            if (len(self.topology.neighbors(node)) == 0):
                self.topology.remove_node(node)

        self.topology = nx.convert_node_labels_to_integers(self.topology)
        return self.topology

    def annotate_edges(self):

        node_ifaces = defaultdict(defaultdict)
        added_edges = set([])

        for edge in self.topology.edges():

            head_id = edge[0]
            tail_id = edge[1]

            if (edge not in added_edges):

                # initialize the node_ifaces dict (if not already)
                if 'prev' not in node_ifaces[head_id]:
                    node_ifaces[head_id]['prev'] = 0
                    node_ifaces[head_id]['ifaces'] = defaultdict(int)

                if 'prev' not in node_ifaces[tail_id]:
                    node_ifaces[tail_id]['prev'] = 0
                    node_ifaces[tail_id]['ifaces'] = defaultdict(int)

                node_ifaces[head_id]['prev'] += 1
                node_ifaces[head_id]['ifaces'][tail_id] = node_ifaces[head_id]['prev']
                # print("added iface %d @ %d for link (%d, %d)" % (node_ifaces[node]['prev'], node, node, link))

                node_ifaces[tail_id]['prev'] += 1
                node_ifaces[tail_id]['ifaces'][head_id] = node_ifaces[tail_id]['prev']
                # print("added iface %d @ %d for link (%d, %d)" % (node_ifaces[link]['prev'], link, link, node))

                added_edges.add(edge)

            self.topology.add_edge(head_id, tail_id,
                e1 = ("%d:%d"   % (head_id, node_ifaces[head_id]['ifaces'][tail_id])), 
                e2 = ("%d:%d"   % (tail_id, node_ifaces[tail_id]['ifaces'][head_id])))

def parse_pop_level_map(rocketfuel_file):

    topology = nx.Graph()
    comment_char = '#'

    # node_ifaces = defaultdict(defaultdict)

    # extract the node names and edges
    node_id = 0
    node_ids = defaultdict(int)
    for line in open(rocketfuel_file, "r").readlines():

        if comment_char in line:
            # split on comment char, keep only the part before
            line, _ = line.split(comment_char, 1)
            line = line.strip()

        if len(line) == 0:
            continue

        try:

            # extract head and tail nodes of the edge
            head = line.split(" -> ")[0].split(",")[0]
            tail = line.split(" -> ")[1].split(",")[0]
            # print("%s -> %s" % (head, tail))

            # add nodes to pop-level topology graph
            if head not in node_ids:
                # add node ids to lookup table
                node_ids[head] = node_id
                node_id += 1
                topology.add_node(node_ids[head], node_str = head)
                # print("added node w/ id %d (str : %s)" % (node_ids[head], head))

            if tail not in node_ids:
                # add node ids to lookup table
                node_ids[tail] = node_id
                node_id += 1
                topology.add_node(node_ids[tail], node_str = tail)
                # print("added node w/ id %d (str : %s)" % (node_ids[tail], tail))

            head_id = node_ids[head]
            tail_id = node_ids[tail]

            topology.add_edge(head_id, tail_id)
            # print("added edge (%d -> %d) (str : %s)" % (head_id, tail_id, str(topology.get_edge_data(head_id, tail_id))))

        except IndexError:
            raise ValueError('Invalid input file. Parsing failed '\
                             'while trying to parse an internal node')

    # the famed Topology object...
    topology_obj = Topology(topology)
    
    # # trim nodes w/ outdegrees > 20 in order to make the evaluation less 
    # # time consuming
    # topology_obj.trim_topology(18)

    # annotate edges w/ iface numbers
    topology_obj.annotate_edges()
    # 
    topology_obj.get_shortest_paths(force = True)
    topology_obj.get_pop_level_statistics(force = True)

    return topology_obj
