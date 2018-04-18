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
from convert_rocketfuel import *

from collections import defaultdict
from collections import OrderedDict

class tp_tree(object):

    def __init__(self):

        self.children = []

        self.router = None
        self.size = None
        self.iface = None

        # FIXME : added for loop prevention when
        # parsing .scn files
        self.riface = None

def get_tp_path_rec(curr_tp_node, tp_path, next_size, size = 1):

    # get sizes of curr_tp_node's children
    sizes = set([])
    for child in curr_tp_node.children:
        sizes.add(child.size)

    # if size is present in sizes, only 
    # follow the children w/ this size from now on
    if size in sizes:
        next_size = size
    else:

        if next_size is None:
            next_size = max(sizes)
        else:
            next_size = next_size

    for child in curr_tp_node.children:

        if child.size == next_size:
            tp_path.append(child.router)
            get_tp_path_rec(child, tp_path, next_size, size)
            break

def get_tp_path(tp_root, size = 1):

    tp_path = [ tp_root.router ]
    get_tp_path_rec(tp_root, tp_path, None, size)
    return tp_path

def extract_tp_path_rec(topology, curr_tp_node, in_iface):

    # find next router
    router = None
    for r in topology.findall('router'):
        if (int(r.get('id')) == curr_tp_node.router):
            router = r
            break

    # get all tps of the router
    tp_ifaces = []
    new_nodes = defaultdict()

    for tps in router.findall('tp'):

        tp_iface    = int(tps.get('iface'))
        tp_size     = int(tps.text)

        # print("""router : %s, tp iface : %s, tp size : %s""" % (router.get('id'), tp_iface, tp_size))
        tp_ifaces.append(tp_iface)

        # for each tp, add a new tp node
        new_tp_node = tp_tree()
        for link in router.findall('link'):

            if int(link.get('local')) == tp_iface:

                new_tp_node.router  = int(link.get('rrouter'))
                new_tp_node.size    = int(tp_size)
                new_tp_node.iface   = int(tp_iface)
                # FIXME : for loop prevention
                new_tp_node.riface  = int(link.get('remote'))

                new_nodes[tp_iface] = new_tp_node
                break

    # for each tp iface, call extract_tp_path_rec() recursively
    for tp_iface in sorted(tp_ifaces):

        # stop recursive calls if next_tp_node has 
        # a local tp entry, i.e., if tp_iface == 0
        if (tp_iface == 0) or (tp_iface == in_iface):
            continue

        # print("added to %s : (%s, %s, %s)" % (curr_tp_node.router, new_nodes[tp_iface].router, new_nodes[tp_iface].iface, new_nodes[tp_iface].size))
        curr_tp_node.children.append(new_nodes[tp_iface])
        extract_tp_path_rec(topology, new_nodes[tp_iface], new_nodes[tp_iface].riface)

def extract_tp_path(scn_file, initial_router):

    topology = et.ElementTree(file = scn_file)
    topology = topology.getroot()

    tp_root = tp_tree()
    tp_root.router = initial_router

    extract_tp_path_rec(topology, tp_root, 0)

    return tp_root

def verify(test_dir):

    for test_file in sorted(glob.glob(os.path.join(test_dir, '*.test'))):
    
        # FIXME: the <test_run> element is useless
        test_run = et.parse(test_file)
        # verify each test within the test run
        for test in test_run.findall('test'):

#             print("""%s: [INFO] test params:\n
# \tid : %s
# \tbf size : %s
# \treq. size : %s
# \tMM mode : %s
# \tIDR mode : %s""" % (
#             sys.argv[0], 
#             test.get('id'), 
#             test.find('bf-size').text, 
#             test.find('req-size').text, 
#             test.find('modes').text.split("-", 1)[0], 
#             test.find('modes').text.split("-", 1)[1]))

            # extract the true and annc. tp paths
            for path in test.find('paths').findall('path'):

                tp_paths = {'true' : [], 'anncd' : []}

                # as announced on .test file
                tp_paths['anncd'] = [ int(x) for x in path.text.split(",") ]
                # as extracted on .scn file
                tp_tree = extract_tp_path(path.get('file'), tp_paths['anncd'][0])
                tp_paths['true'] = get_tp_path(tp_tree)

                if tp_paths['true'] != tp_paths['anncd']:
                    print("differences in tp paths [ size %d: ]" % (1))
                    print("\tfile: %s"  % (path.get('file')))                    
                    print("\tanncd: %s" % (tp_paths['anncd']))
                    print("\ttrue: %s"  % (tp_paths['true']))

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--print-stats", 
         help = """print statistics about RocketFuel dataset""",
         action = "store_true")

    parser.add_argument(
        "--test-dir", 
         help = """dir on which to save .test and .scn files""")

    args = parser.parse_args()

    if not args.test_dir:
        sys.stderr.write("""%s: [ERROR] please supply a .test dir\n""" % sys.argv[0])
        parser.print_help()
        sys.exit(1)

    verify(args.test_dir)