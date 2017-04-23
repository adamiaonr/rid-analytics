import argparse
import re
import sys
import glob
import os
import math
import binascii

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import xml.etree.cElementTree as et

from xml.dom import minidom
from prettytable import PrettyTable

from collections import defaultdict
from collections import OrderedDict

# post simple statistics about each AS-level topology:
#   - nr. of nodes
#   - nr. of links
#   - median outdegree
def print_statistics(rocketfuel_dir):

    stats = OrderedDict()

    for filename in sorted(glob.glob(os.path.join(rocketfuel_dir, '*r0.cch'))):

        topology_id = filename.split("/")[-1]
        topology_id = topology_id.split(".", 1)[0]
        topology = parse_isp_map(filename)

        stats[topology_id] = defaultdict()

        stats[topology_id]['nodes'] = topology.number_of_nodes()
        stats[topology_id]['links'] = topology.number_of_edges()
        degree_sequence = sorted(nx.degree(topology).values(), reverse = True)
        stats[topology_id]['outdegree'] = np.median(degree_sequence)

    table = PrettyTable(['as id', '# nodes', '# links', 'median outdegree'])

    for topology_id in stats:
        table.add_row([
            topology_id,
            stats[topology_id]['nodes'],
            stats[topology_id]['links'],
            stats[topology_id]['outdegree']
            ])

    print("")
    print(table)
    print("")

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
            topology = parse_pop_level_map(filename)

            stats[topology_id] = defaultdict()

            stats[topology_id]['nodes'] = topology.number_of_nodes()
            stats[topology_id]['links'] = topology.number_of_edges()
            degree_sequence = sorted(nx.degree(topology).values(), reverse = True)
            print(degree_sequence)
            stats[topology_id]['min-outdegree'] = np.min(degree_sequence)
            stats[topology_id]['med-outdegree'] = np.median(degree_sequence)
            stats[topology_id]['max-outdegree'] = np.max(degree_sequence)

    table = PrettyTable(['as id', '# nodes', '# links', 'min. outdegree', 'median outdegree', 'max. outdegree'])

    for topology_id in stats:
        table.add_row([
            topology_id,
            stats[topology_id]['nodes'],
            stats[topology_id]['links'],
            stats[topology_id]['min-outdegree'],
            stats[topology_id]['med-outdegree'],
            stats[topology_id]['max-outdegree']
            ])

    print("")
    print(table)
    print("")

# add true positive entries to a topology (.scn file). we start by applying 
# changes in a networkx object, then call convert_to_scn() to generate a .scn 
# file.
#
# algorithm:
#   - find all shortest paths in topology originated on content_source_id 
#   - add entries of size tp_size on links in shortest path
def add_content_route(topology, content_source_id, tp_size = 1, radius = -1):

    routes = nx.single_source_shortest_path(topology, content_source_id)

    # add a tp entry to the local iface of content_source_id. this is basically 
    # an edge with the same head and tail.
    topology.add_edge(content_source_id, content_source_id, 
        type = 'internal', 
        e1 = ("%d:%d"   % (content_source_id, 0)), 
        e2 = ("%d:%d"   % (content_source_id, 0)),
        tp = ("%d:%d:%d" % (content_source_id, 0, tp_size)))

    # iterate through the links in the path, so that we add tp information 
    # to them
    for route in routes:
        for i in xrange(len(routes[route]) - 1):
            
            # get src and dst ends of the edge
            src, dst = routes[route][i], routes[route][i + 1]
            # determine the interface of dst to which we should add 
            # the tp info
            edge = topology.get_edge_data(src, dst);
            iface = 0
            if int(edge['e1'].split(':', 1)[0]) == dst:
                iface = edge['e1']
            else:
                iface = edge['e2']

            # the tp info is a tuple of the form <router>:<outgoing iface>:<size> 
            topology.add_edge(src, dst, tp = ("%s:%d" % (iface, tp_size)))

    return 0

def get_shortest_paths(topology):

    # compute shortest paths from every node to every other node. the paths 
    # are stored as strings, sequences of node ids, separated by ','. the structure 
    # of the returned dictionary is:
    #
    #   - [source 1] -> [ 
    #                   "<source 1>,<path node 1>,<path node 2>, ... ,<node 0>",
    #                   "<source 1>,<path node 1>,<path node 2>, ... ,<node 1>",
    #                   (...),
    #                   "<source 1>,<path node 1>, ... ,<node n>"
    #                   ]
    #
    #   - [source 2] -> [ 
    #                   "<source 2>,<path node 1>,<path node 2>, ... ,<node 0>",
    #                   "<source 2>,<path node 1>,<path node 2>, ... ,<node 1>",
    #                   (...),
    #                   "<source 2>,<path node 1>, ... ,<node n>"
    #                   ]
    #   - (...)
    shortest_paths = defaultdict(list)
    for source in topology.nodes():

        routes = nx.single_source_shortest_path(topology, source)

        for dst, route in routes.iteritems():
            # we transform the list of nodes in the route into a string, which 
            # then allows for more convenient search on get_iface_distr()
            route = ','.join(("%03d" % (node)) for node in route)
            shortest_paths[source].append(route)

    return shortest_paths

# FIXME : this algorithm is terribly inefficient...
def get_fwd_dist(topology, shortest_paths, router):

    # fwd_dist values
    fwd_dist = defaultdict(float)
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
        
        # if (i + 1) == neighbor_num:
        #     fwd_dist[neighbor] = float(node_num - 1) - float(len(visited_sources))
        #     continue
        search_str = ("%03d,%03d" % (neighbor, router))
        
        for source in shortest_paths:

            if source in visited_sources:
                continue;

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

def to_tree_bitmask(topology, iface_trees):

    # calculate nr. of bytes of the bitmask. note we can represent 8 nodes per 
    # byte
    nr_bytes = int(math.ceil(float(len(topology.nodes())) / 8.0))
    # initialize tree bitmask to 0s
    tree_bitmask = [0] * nr_bytes

    for source in iface_trees:
        byte_index = int((source / 8))
        # print("source : %d, tree_bitmask[%d] = %d" % (source, byte_index, tree_bitmask[byte_index]))
        tree_bitmask[byte_index] |= (1 << (source % 8))
        # print("bit index = %d, tree_bitmask[%d] = %d" % ((source % 8), byte_index, tree_bitmask[byte_index]))

    return nr_bytes, binascii.hexlify(bytearray(tree_bitmask))

def convert_to_scn(topology, req_size = 15):

    shortest_paths = get_shortest_paths(topology)
    topology_block = et.Element("topology")

    # cycle through each node in the topology. we then write 5 diff. types of 
    # blocks: 
    #   - <router id=""> block
    #   - <link local="" remote="" rrouter=""> block
    #   - <fwd_size_dist size="">x.xx</fwd_size_dist> block
    #   - <fwd_dist iface="y">x.xx</fwd_dist>
    #   - <tp iface="y">x</tp>
    for router in topology.nodes():

        added_local_iface = False

        # add <router id=""> block
        router_block = et.SubElement(topology_block, "router", id = str(router))

        # get the fwd_dist values for the router's ifaces
        fwd_dist, iface_trees = get_fwd_dist(topology, shortest_paths, router)

        # add <link> blocks, one for each neighbor of router
        for i, neighbor_router in enumerate(topology.neighbors(router)):

            # extract the attributes of the link
            edge = topology.get_edge_data(router, neighbor_router);

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

            # add <fwd_size_dist> blocks, indicating the distribution of 
            # forwarding entry sizes at each link.
            for i in xrange(req_size):
                fwd_size_dist_block = et.SubElement(link_block, "fwd_size_dist", size = str(i + 1)).text = "0.00"

            # add <tp> block (if edge has 'tp' attribute)
            # <tp iface="2">1</tp>
            if ('tp' in edge) and (int(edge['tp'].split(':')[0]) == router):
                tp_block = et.SubElement(router_block, "tp", iface = edge['tp'].split(':')[1]).text = edge['tp'].split(':')[2]

            # add tree_bitmask
            nr_bytes, tree_bitmask = to_tree_bitmask(topology, list(iface_trees[neighbor_router]))
            tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", size = str(nr_bytes)).text = tree_bitmask

            if (neighbor_router == router):
                added_local_iface = True

        # add a special local link (represents local delivery), if not added 
        # already
        if not added_local_iface:
            link_block = et.SubElement(router_block, "link", 
                local   = str(0), 
                remote  = str(0), 
                rrouter = str(router))

            for i in xrange(req_size):
                fwd_size_dist_block = et.SubElement(link_block, "fwd_size_dist", size = str(i + 1)).text = "0.00"

            nr_bytes, tree_bitmask = to_tree_bitmask(topology, [router])
            tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", size = str(nr_bytes)).text = tree_bitmask

        # add <fwd_dist> blocks for each iface
        # if the local iface isn't calculated in get_fwd_dist() (e.g. that 
        # only happens if the link (router, router) has been added), add it 
        # now.
        if not added_local_iface:
            fwd_dist_block = et.SubElement(router_block, "fwd_dist", iface = str(0)).text = ("%.4f" % (fwd_dist[router]))
        # finally, the ifaces pointing to each neighboring router
        for neighbor_router in topology.neighbors(router):

            # find iface which links to neighbor_router
            edge = topology.get_edge_data(router, neighbor_router)
            local_iface = 0
            if int(edge['e1'].split(':', 1)[0]) == router:
                local_iface = edge['e1'].split(':', 1)[1]
            else:
                local_iface = edge['e2'].split(':', 1)[1]

            fwd_dist_block = et.SubElement(router_block, "fwd_dist", iface = str(local_iface)).text = ("%.4f" % (fwd_dist[neighbor_router]))


    xmlstr = minidom.parseString(et.tostring(topology_block)).toprettyxml(indent="    ")
    with open("topology.scn", "w") as f:
        f.write(xmlstr)

def draw_pop_level_map(topology):

    pos = nx.spring_layout(topology)
    nx.draw_networkx_nodes(
        topology, 
        pos,
        node_color = 'red',
        node_size = 500,
        alpha = 0.5)

    nx.draw_networkx_edges(
        topology,
        pos,
        width = 1.0, 
        alpha = 0.5)

    labels={}
    for router in topology.nodes():
        labels[router] = str(router)
    nx.draw_networkx_labels(topology, pos, labels, font_size = 10)

    # save the figure in <rocketfuel-file>.pdf
    plt.savefig(args.data_path.split("/")[-2] + ".pdf", bbox_inches='tight', format = 'pdf')

def parse_pop_level_map(rocketfuel_file):

    topology = nx.Graph()
    comment_char = '#'

    node_ifaces = defaultdict(defaultdict)

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

            if (topology.has_edge(head_id, tail_id) == False):

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

            topology.add_edge(head_id, tail_id,
                e1 = ("%d:%d"   % (head_id, node_ifaces[head_id]['ifaces'][tail_id])), 
                e2 = ("%d:%d"   % (tail_id, node_ifaces[tail_id]['ifaces'][head_id])))

            # print("added edge (%d -> %d) (str : %s)" % (head_id, tail_id, str(topology.get_edge_data(head_id, tail_id))))

        except IndexError:
            raise ValueError('Invalid input file. Parsing failed '\
                             'while trying to parse an internal node')

    return topology

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--data-path", 
         help = """path w/ topology files (e.g. .cch for isp maps, edge folders 
            for pop-level maps, etc.)""")

    parser.add_argument(
        "--add-tp-source", 
         help = """e.g. '--add-tp-source <source id>:<size>:<radius>'""")

    parser.add_argument(
        "--print-stats", 
         help = """print statistics about RocketFuel dataset""",
         action = "store_true")

    args = parser.parse_args()

    # quit if a dir w/ causality files hasn't been provided
    if not args.data_path:
        sys.stderr.write("""%s: [ERROR] please supply a data path\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    if args.print_stats:
        print_pop_level_statistics(args.data_path)
        # print_statistics(args.data_path)
        sys.exit(0)

    # build a networkx topology out of a Rocketfuel pop-level topology
    topology = parse_pop_level_map(args.data_path)    
    draw_pop_level_map(topology)

    # if requested, add a true positive content source
    if args.add_tp_source:
        add_content_route(
            topology, 
            int(args.add_tp_source.split(':')[0]), 
            int(args.add_tp_source.split(':')[1]))

    # finally, convert the topology to an .scn file
    convert_to_scn(topology)
