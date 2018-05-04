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

def verify_paths(test_dir):

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

def verify_tps(test_dir):

    for test_file in sorted(glob.glob(os.path.join(test_dir, '*.test'))):
    
        # FIXME: the <test_run> element is useless
        test_run = et.parse(test_file)
        # verify each test within the test run
        for test in test_run.findall('test'):

            infractions = defaultdict(list)
            for path in test.find('paths').findall('path'):

                p = [ int(x) for x in path.text.split(",") ]

                scn_file = path.get('file')
                topology = et.ElementTree(file = scn_file)
                topology = topology.getroot()

                for r in topology.findall('router'):

                    tps = defaultdict(list)
                    for tp in r.findall('tp'):
                        tps[int(tp.text)].append(tp)

                    for tp_size in tps:
                        if len(tps[tp_size]) != 1:
                            infractions[("%s-%d-%d:" % (test.get('id'), p[0], p[-1]))].append(("\t%d : %d tps of size %d" % (int(r.get('id')), len(tps[tp_size]), tp_size)))

            for infraction in infractions:
                print("%s" % (infraction))
                for s in infractions[infraction]:
                    print("%s" % (s))

# def remove_from_bitmask(router, iface, size, dst):

#     byte_pos = int(dst / 8)
#     bit_pos  = int(dst % 8)

#     for l in router.findall('link'): 

#         if int(l.get('local')) != iface:
#             continue

#         for tb in l.findall('tree_bitmask'):

#             if int(tb.get('size')) == size:


def fix_tps(test_dir, size = 1):

    for test_file in sorted(glob.glob(os.path.join(test_dir, '*.test'))):
    
        # FIXME: the <test_run> element is useless
        test_run = et.parse(test_file)
        # verify each test within the test run
        for test in test_run.findall('test'):

            infractions = defaultdict(list)
            for path in test.find('paths').findall('path'):

                scn_file = path.get('file')
                topology = et.ElementTree(file = scn_file)
                topology = topology.getroot()

                altered = False

                for r in topology.findall('router'):

                    # if # of tps > 1, remove the 'extra' tp
                    tps = len(r.findall('tp'))
                    if tps > 1:

                        _path = [ int(x) for x in path.text.split(",") ]

                        print("@test : %s" % (scn_file.split("/")[-1].rstrip('.scn')))
                        print("@router : %d" % (int(r.get('id'))))

                        # if the curr router is in the tp path from 
                        # src to dst, remove the 'extra' tp accordingly
                        if int(r.get('id')) in _path:

                            # find which is the 'correct' next router
                            r_pos = [i for i, x in enumerate(_path) if x == int(r.get('id'))][0]
                            next_router = _path[r_pos + 1]
                            print("%d -> %d in %s" % (_path[r_pos], _path[r_pos + 1], _path))

                            # find the iface which points to next_router
                            next_router_iface = 0
                            for l in r.findall('link'):

                                if int(l.get('rrouter')) == next_router:
                                    next_router_iface = int(l.get('local'))
                                    break

                            print("%d[%d] -> %d in %s" % (_path[r_pos], next_router_iface, _path[r_pos + 1], _path))

                            # remove the tps which don't point towards 
                            for tp in r.findall('tp'):

                                if int(tp.text) != size:
                                    continue

                                if int(tp.get('iface')) != next_router_iface:
                                    print("will remove %d:%d:%d (vs. %d:%d:%d)" % (int(r.get('id')), int(tp.get('iface')), int(tp.text), int(r.get('id')), next_router_iface, size))
                                    r.remove(tp)
                                    # remove_from_bitmask(r, iface, int(tp.get('iface')))
                                    altered = True

                        # else, just remove all but one of them
                        else:

                            tps = defaultdict(list)
                            for tp in r.findall('tp'):
                                tps[int(tp.text)].append(tp)

                            if size in tps:

                                tp_num = len(tps[size])
                                for tp in tps[size]:

                                    if tp_num > 1:
                                        print("will remove %d:%d" % (int(r.get('id')), int(tp.get('iface'))))
                                        r.remove(tp)
                                        tp_num = tp_num - 1
                                        altered = True
                                    
                                    else:
                                        break

                # finally, overwrite the .scn file
                if altered:
                    xmlstr = minidom.parseString(et.tostring(topology).replace('    ', '').replace('\t', '').replace('\r', '').replace('\n', '')).toprettyxml()
                    with open(scn_file, "w") as f:
                        f.write(xmlstr)

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

    parser.add_argument(
        "--case", 
         help = """verification or 'fix' case""")

    parser.add_argument(
        "--size", 
         help = """size of tp verification or fix (default: 1)""")

    args = parser.parse_args()

    if not args.test_dir:
        sys.stderr.write("""%s: [ERROR] please supply a .test dir\n""" % sys.argv[0])
        parser.print_help()
        sys.exit(1)

    if args.case == 'verify-tps':
        verify_tps(args.test_dir)

    elif args.case == 'fix-tps':

        size = 1
        if args.size:
            size = int(args.size)

        fix_tps(args.test_dir, size = size)
        verify_tps(args.test_dir)

    elif args.case == 'verify-paths':
        verify_paths(args.test_dir)

    else:
        sys.stderr.write("""%s: [ERROR] please supply a valid case""" % (sys.argv[0]))
        parser.print_help()
        sys.exit(1)
