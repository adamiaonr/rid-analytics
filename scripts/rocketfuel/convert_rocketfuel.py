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

DEFAULT_DIFF_DISTR_DIR = "/home/adamiaonr/workbench/rid-analytics/experiments/examples/diff-distributions"

def get_pop_level_statistics(topology):

        stats = defaultdict()

        stats['nodes'] = topology.number_of_nodes()
        stats['links'] = topology.number_of_edges()
        stats['outdegree-list'] = nx.degree(topology)
        degree_sequence = sorted(nx.degree(topology).values(), reverse = True)

        avg_path_outdegrees = []
        path_lengths = []
        shortest_paths = get_shortest_paths(topology)
        for source in shortest_paths:
            for path in shortest_paths[source]:

                path_lengths.append(len(path.split(",")))
                path_routers = [int(r) for r in path.split(",")]
                avg_path_outdegrees.append(float(sum([stats['outdegree-list'][r] for r in path_routers])) / (float(len(path.split(",")))))

        stats['avg-path-outdegree'] = list(set(avg_path_outdegrees))
        stats['path-lengths'] = sorted(list(set(path_lengths)), reverse = True)

        stats['min-outdegree'] = np.min(degree_sequence)
        stats['med-outdegree'] = np.median(degree_sequence)
        stats['max-outdegree'] = np.max(degree_sequence)

        return stats

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

            stats[topology_id] = get_pop_level_statistics(topology)

    table = PrettyTable(['as id', '# nodes', '# links', 'min. outdegree', 'median outdegree', 'max. outdegree', 'diameter'])

    for topology_id in stats:
        table.add_row([
            topology_id,
            stats[topology_id]['nodes'],
            stats[topology_id]['links'],
            stats[topology_id]['min-outdegree'],
            stats[topology_id]['med-outdegree'],
            stats[topology_id]['max-outdegree'],
            (stats[topology_id]['path-lengths'][0] - 1)
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

def add_example_path(example_paths, length_cat, outdegree_cat, path_routers):

    if length_cat not in example_paths:
        example_paths[length_cat] = defaultdict(list)

    example_paths[length_cat][outdegree_cat].append(path_routers)

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
        # print(shortest_paths[source])

    return shortest_paths

def generate_path_examples(topology, pick = (-1, -1)):

    # sample a small set of path examples to use in experiments. we vary the 
    # topology aspects which affect efficiency and correctness:
    #
    #   1) nr. of hops, or path length : the more hops, the higher the chance of 
    #                                    forwarding errors
    #   2) outdegree    : the larger the outdegree, the larger the number of 
    #                     links over which fps matches are forwarded
    #
    # as such, we choose paths w/ the following combinations:
    #   1 2 3) short, low / median / high outdegree
    #   4 5 6) median length, low / median / high outdegree
    #   7 8 9) long length, low / median / high outdegree
    #
    # choosing 2 per example, we have a max. of 18 representative paths, which 
    # should be ok for testing time purposes on a topology Telstra or AT&T
    path_examples = defaultdict()
    # get statistics of topology
    stats = get_pop_level_statistics(topology)
    shortest_paths = get_shortest_paths(topology)
    print("outdegrees : %s" % (sorted(set(nx.degree(topology).values()), reverse = True)))
    # # get outdegree sums of topologies
    # print("path outdegree sums : %s" % (sorted(stats['outdegree-sums'], reverse = True)))

    avg_path_outdegrees = defaultdict()
    avg_path_outdegrees['short'] = defaultdict(list)
    avg_path_outdegrees['median'] = defaultdict(list)
    avg_path_outdegrees['long'] = defaultdict(list)


    for source in reversed([s for s in shortest_paths]):
        for path in shortest_paths[source]:

            path_length = len(path.split(","))
            # don't gather paths for lengths equal to 1
            if path_length < 2:
                continue
            path_routers = [int(r) for r in path.split(",")]
            avg_path_outdegree = float(sum([stats['outdegree-list'][r] for r in path_routers])) / float(path_length)

            # short lengths
            if path_length < np.percentile(stats['path-lengths'], 25):

                if avg_path_outdegree < np.percentile(stats['avg-path-outdegree'], 25):
                    add_example_path(path_examples, 'short', 'low', path_routers)
                    avg_path_outdegrees['short']['low'].append(avg_path_outdegree)
                elif avg_path_outdegree > np.percentile(stats['avg-path-outdegree'], 40) and avg_path_outdegree < np.percentile(stats['avg-path-outdegree'], 60):
                    add_example_path(path_examples, 'short', 'median', path_routers)
                    avg_path_outdegrees['short']['median'].append(avg_path_outdegree)
                elif avg_path_outdegree > np.percentile(stats['avg-path-outdegree'], 75):
                    add_example_path(path_examples, 'short', 'high', path_routers)
                    avg_path_outdegrees['short']['high'].append(avg_path_outdegree)

            # median lengths
            if path_length > np.percentile(stats['path-lengths'], 40) and path_length < np.percentile(stats['path-lengths'], 60):

                if avg_path_outdegree < np.percentile(stats['avg-path-outdegree'], 25):
                    add_example_path(path_examples, 'median', 'low', path_routers)
                    avg_path_outdegrees['median']['low'].append(avg_path_outdegree)
                elif avg_path_outdegree > np.percentile(stats['avg-path-outdegree'], 40) and avg_path_outdegree < np.percentile(stats['avg-path-outdegree'], 60):
                    add_example_path(path_examples, 'median', 'median', path_routers)
                    avg_path_outdegrees['median']['median'].append(avg_path_outdegree)
                elif avg_path_outdegree > np.percentile(stats['avg-path-outdegree'], 75):
                    add_example_path(path_examples, 'median', 'high', path_routers)
                    avg_path_outdegrees['median']['high'].append(avg_path_outdegree)

            # long lengths
            if path_length > np.percentile(stats['path-lengths'], 75):

                if avg_path_outdegree < np.percentile(stats['avg-path-outdegree'], 25):
                    add_example_path(path_examples, 'long', 'low', path_routers)
                    avg_path_outdegrees['long']['low'].append(avg_path_outdegree)
                elif avg_path_outdegree > np.percentile(stats['avg-path-outdegree'], 40) and avg_path_outdegree < np.percentile(stats['avg-path-outdegree'], 60):
                    add_example_path(path_examples, 'long', 'median', path_routers)
                    avg_path_outdegrees['long']['median'].append(avg_path_outdegree)
                elif avg_path_outdegree > np.percentile(stats['avg-path-outdegree'], 75):
                    add_example_path(path_examples, 'long', 'high', path_routers)
                    avg_path_outdegrees['long']['high'].append(avg_path_outdegree)


    # chosen_path = defaultdict()
    # chosen_path_avg_outdegree = defaultdict()

    # for path_length in path_examples:

    #     chosen_path[path_length] = defaultdict()
    #     chosen_path_avg_outdegree[path_length] = defaultdict()

    #     for path_outdegree in path_examples[path_length]:

    #         # nr_paths = len(path_examples[path_length][path_outdegree])

    #         if pick != (-1, -1):
    #             for i, path in enumerate(path_examples[path_length][path_outdegree]):
    #                 if path[0] == pick[0] and path[-1] == pick[1]:
    #                     chosen_path[path_length][path_outdegree] = path
    #                     chosen_path_avg_outdegree[path_length][path_outdegree] = avg_path_outdegrees[path_length][path_outdegree][i]

    #         else:
    #             chosen_path[path_length][path_outdegree] = path_examples[path_length][path_outdegree][0]
    #             chosen_path_avg_outdegree[path_length][path_outdegree] = avg_path_outdegrees[path_length][path_outdegree][0]

    # return chosen_path, chosen_path_avg_outdegree
    return path_examples, avg_path_outdegrees

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
    # print("got %d bytes" % (nr_bytes))
    # initialize tree bitmask to 0s
    tree_bitmask = [0] * nr_bytes

    for source in iface_trees:
        byte_index = int((source / 8))
        # print("source %d, byte index : %d" % (source, byte_index))
        # print("source : %d, tree_bitmask[%d] = %d" % (source, byte_index, tree_bitmask[byte_index]))
        tree_bitmask[byte_index] |= (1 << (source % 8))
        # print("bit index = %d, tree_bitmask[%d] = %d" % ((source % 8), byte_index, tree_bitmask[byte_index]))

    return nr_bytes, binascii.hexlify(bytearray(tree_bitmask))

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

def convert_to_scn(topology, entry_sizes, table_size = 100000000, req_size = 5, scn_filename = 'topology.scn'):

    topology_block = et.Element("topology")

    # get shortest path list and determine a ttl value for the 
    # topology (2 * longest path, as in http://www.map.meteoswiss.ch/map-doc/ftp-probleme.htm)
    shortest_paths = get_shortest_paths(topology)

    ttl = 0
    for source in shortest_paths:
        for path in shortest_paths[source]:

            path_length = (len(path.split(",")) - 1)
            if ttl < path_length:
                ttl = path_length

    ttl_block = et.SubElement(topology_block, "ttl").text = str(ttl)

    # load |F\R| distributions
    load_request_entry_diff_distributions(req_size, topology_block)
    # fwd table size
    fwd_table_size_block = et.SubElement(topology_block, "fwd_table_size").text = str(table_size)

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

                fwd_size_dist_block_text = "0.00"
                if str(i + 1) in entry_sizes:
                    fwd_size_dist_block_text = ("%.2f" % (entry_sizes[str(i + 1)]))

                fwd_size_dist_block = et.SubElement(link_block, "fwd_size_dist", size = str(i + 1)).text = fwd_size_dist_block_text

            # add <tp> block (if edge has 'tp' attribute)
            # <tp iface="2">1</tp>
            if ('tp' in edge) and (int(edge['tp'].split(':')[0]) == router):
                tp_block = et.SubElement(router_block, "tp", iface = edge['tp'].split(':')[1]).text = edge['tp'].split(':')[2]

            # add tree_bitmask
            if (neighbor_router == router):
                nr_bytes, tree_bitmask = to_tree_bitmask(topology, [router])
                tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", size = str(nr_bytes)).text = tree_bitmask

                added_local_iface = True
            
            else:
                nr_bytes, tree_bitmask = to_tree_bitmask(topology, list(iface_trees[neighbor_router]))
                tree_bitmask_block = et.SubElement(link_block, "tree_bitmask", size = str(nr_bytes)).text = tree_bitmask


        # add a special local link (represents local delivery), if not added 
        # already
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
    with open(scn_filename, "w") as f:
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

def trim_topology(topology, max_outdegree):

    # extract stats from current topology
    stats = get_pop_level_statistics(topology)

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

            for i, neighbor_router in enumerate(topology.neighbors(router)):

                if to_trim == 0:
                    break

                if stats['outdegree-list'][neighbor_router] == r:
                    topology.remove_node(neighbor_router)
                    to_trim -= 1

    # remove nodes with no neighbors
    for node in topology.nodes():

        if (len(topology.neighbors(node)) == 0):
            topology.remove_node(node)

    topology = nx.convert_node_labels_to_integers(topology)
    return topology

def annotate_edges(topology):

    node_ifaces = defaultdict(defaultdict)
    added_edges = set([])

    for edge in topology.edges():

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

        topology.add_edge(head_id, tail_id,
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

    # trim nodes w/ outdegrees > 20 in order to make the evaluation less 
    # time consuming
    topology = trim_topology(topology, 18)
    # annotate edges w/ iface numbers
    annotate_edges(topology)

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
        "--req-size", 
         help = """max nr. of prefix components which can be encoded in a BF""")

    parser.add_argument(
        "--add-tp-source", 
         help = """e.g. '--add-tp-source <source id>:<size>:<radius>'""")

    parser.add_argument(
        "--table-size", 
         help = """e.g. '--table-size 10000000'""")

    parser.add_argument(
        "--entry-sizes", 
         help = """e.g. '--entry-sizes <entry-size>:<size %%>|<entry-size>:<size %%>|...|<entry-size>:<size %%>'""")

    parser.add_argument(
        "--print-stats", 
         help = """print statistics about RocketFuel dataset""",
         action = "store_true")

    args = parser.parse_args()

    entry_sizes = defaultdict();
    if args.entry_sizes:
        entry_size_records = args.entry_sizes.split("|")
        for record in entry_size_records:
            entry_sizes[record.split(":", 1)[0]] = (float(record.split(":", 1)[1]) / 100.0)

    req_size = 5
    if args.req_size:
        req_size = int(args.req_size)

    table_size = 10000000
    if args.table_size:
        table_size = int(args.table_size)

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
    generate_path_examples(topology)

    # if requested, add a true positive content source
    if args.add_tp_source:
        add_content_route(
            topology, 
            int(args.add_tp_source.split(':')[0]), 
            int(args.add_tp_source.split(':')[1]))

    # finally, convert the topology to an .scn file
    convert_to_scn(topology, entry_sizes, table_size, req_size)
