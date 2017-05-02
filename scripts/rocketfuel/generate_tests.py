import argparse
import re
import sys
import glob
import os
import math
import binascii
import random
import convert_rocketfuel

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import xml.etree.cElementTree as et

from xml.dom import minidom
from prettytable import PrettyTable

from collections import defaultdict
from collections import OrderedDict

def add_new_test(main_block, topology, test_parameters, path_examples, avg_outdegrees, test_dir):

    # a test run is composed by:
    #   1) a topology nr. (e.g. 4755)
    #   2) bf size (in bit)
    #   3) a request size
    #   4) fwd entry size distribution
    #   5) table size
    #   6) modes for mm, eh, resolution
    #   7) a true positive path (source and destination nodes)

    # format of test id is : <topology nr>:<bf-size>:<request size>:<fwd entry size>:<table size>:<modes>:<source>:<destination>
    test_id = ("%d-%03d-%02d-%02d-%d-%s" % (
        test_parameters['topology-nr'], 
        test_parameters['bf-size'],
        test_parameters['req-size'],
        [int(x) for x in test_parameters['entry-sizes']][0],
        test_parameters['table-size'],
        test_parameters['modes']))

    test_block = et.SubElement(main_block, "test", id = test_id)
    
    for param in test_parameters:

        if param == 'entry-sizes':

            entry_sizes_list = []
            for e in test_parameters['entry-sizes']:
                entry_sizes_list.append(("%02d:%.2f" % (int(e), test_parameters['entry-sizes'][e])))

            param_block = et.SubElement(test_block, param).text = '|'.join(s for s in entry_sizes_list)

        else:
            param_block = et.SubElement(test_block, param).text = str(test_parameters[param])

    paths_block = et.SubElement(test_block, "paths")
    for path_length in path_examples:

        if path_length != 'median':
            continue

        for path_outdegree in path_examples[path_length]:

            if path_outdegree != 'median':
                continue

            # save the scenario info about the particular test
            path = path_examples[path_length][path_outdegree]
            
            # generate the .scn filename for this scenario
            scn_filename = ("configs/%s-%03d-%03d.scn" % (test_id, path[0], path[-1]))
            scn_filename = os.path.join(test_dir, scn_filename)

            new_path_block = et.SubElement(paths_block, "path", length = path_length, outdegre = path_outdegree, avg_outdegree = str(avg_outdegrees[path_length][path_outdegree]), file = scn_filename).text = str(','.join(("%03d" % (x)) for x in path))

            # add true positives
            _topology = topology.copy()
            convert_rocketfuel.add_content_route(
                _topology, 
                path[-1], 
                int([x for x in test_parameters['entry-sizes']][0]))

            # now save the actual .scn file
            convert_rocketfuel.convert_to_scn(
                _topology, 
                test_parameters['entry-sizes'], 
                test_parameters['table-size'], 
                test_parameters['req-size'], 
                scn_filename)

    results_dir_block = et.SubElement(test_block, "results_dir").text = os.path.join(test_dir, "results/")

    return 0

def generate_test(test_dir, topology_file, req_sizes, bf_sizes, entry_size_records, table_sizes, modes, picks):

    main_block = et.Element("test_run")

    # extract the topology nr.
    topology_nr = int(topology_file.split("/")[-2])
    topology_nr_block = et.SubElement(main_block, "topology", file = topology_file).text = ("%d" % (topology_nr))

    # build the basic topology (no .scn file yet), and 
    # generate path examples
    topology = convert_rocketfuel.parse_pop_level_map(topology_file)
    path_examples, avg_outdegrees = convert_rocketfuel.generate_path_examples(topology)

    chosen_path = defaultdict()
    chosen_path_avg_outdegree = defaultdict()

    for path_length in path_examples:

        chosen_path[path_length] = defaultdict()
        chosen_path_avg_outdegree[path_length] = defaultdict()

        for path_outdegree in path_examples[path_length]:

            # nr_paths = len(path_examples[path_length][path_outdegree])
            picks_key = ("%s:%s" % (path_length, path_outdegree))
            if picks_key in picks:
                for i, path in enumerate(path_examples[path_length][path_outdegree]):

                    if (int(path[0]) == picks[picks_key][0]) and (int(path[-1]) == picks[picks_key][1]):
                        chosen_path[path_length][path_outdegree] = path
                        chosen_path_avg_outdegree[path_length][path_outdegree] = avg_outdegrees[path_length][path_outdegree][i]

            else:
                chosen_path[path_length][path_outdegree] = path_examples[path_length][path_outdegree][0]
                chosen_path_avg_outdegree[path_length][path_outdegree] = avg_outdegrees[path_length][path_outdegree][0]

    # build the test cases w/ good ol' nested for loops
    for req_size in req_sizes:
        for bf_size in bf_sizes:
            for record in entry_size_records:

                entry_sizes = defaultdict();
                entry_size_proportions = record.split("|")

                for p in entry_size_proportions:
                    entry_sizes[p.split(":", 1)[0]] = (float(p.split(":", 1)[1]) / 100.0)
                print("entry-sizes = %s" % str(entry_sizes))

                for table_size in table_sizes:
                    for mode in modes:

                        test_parameters = defaultdict()

                        test_parameters['topology-nr'] = topology_nr
                        test_parameters['req-size'] = req_size
                        test_parameters['bf-size'] = bf_size
                        test_parameters['entry-sizes'] = entry_sizes
                        test_parameters['table-size'] = int(table_size)
                        test_parameters['modes'] = mode.replace(":", "-")

                        add_new_test(main_block, topology, test_parameters, chosen_path, chosen_path_avg_outdegree, test_dir)

    xmlstr = minidom.parseString(et.tostring(main_block)).toprettyxml(indent="    ")
    with open(os.path.join(test_dir, ("configs/%d.test" % (topology_nr))), "w") as f:
        f.write(xmlstr)

    return 0

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--test-dir", 
         help = """dir on which to save .test and .scn files""")

    parser.add_argument(
        "--topology-file", 
         help = """filename of edge.wt topology file""")

    parser.add_argument(
        "--req-sizes", 
         help = """max nr. of prefix components which can be encoded in a BF. 
            e.g. '--req-sizes 5:10:15'""")

    parser.add_argument(
        "--bf-sizes", 
         help = """size of BFs in bit. e.g. '--bf-sizes 192:256'""")

    parser.add_argument(
        "--entry-sizes", 
         help = """e.g. '--entry-sizes <entry-size>:<size %%>|<entry-size>:<size %%>,<entry-size>:<size %%>|...'""")

    parser.add_argument(
        "--table-sizes", 
         help = """e.g. '--table-sizes 10000000:100'""")

    parser.add_argument(
        "--modes", 
         help = """MM_MODE (0 for 'flood', 1 for 'random', or 2 for 'fallback'), RES_MODE (0 for 'drop packets', 1 for 'resolve w/ fallback'), in that order, separated by ':' and '|' e.g. '--modes 0:0|0:1'""")

    parser.add_argument(
        "--pick", 
         help = """e.g. '--pick long:median:36:42|short:median:5:78'""")    

    args = parser.parse_args()

    req_sizes = []
    if args.req_sizes:
        req_sizes = [int(x) for x in args.req_sizes.split(":")]

    bf_sizes = []
    if args.bf_sizes:
        bf_sizes = [int(x) for x in args.bf_sizes.split(":")]

    entry_size_records = []
    if args.entry_sizes:
        entry_size_records = args.entry_sizes.split(",")

    table_sizes = []
    if args.table_sizes:
        table_sizes = args.table_sizes.split(":")

    modes = []
    if args.modes:
        modes = args.modes.split("|")

    picks = defaultdict()
    if args.pick:

        pick = args.pick.split("|")
        for p in pick:
            pp = p.split(":")
            picks[("%s:%s" % (pp[0], pp[1]))] = (int(pp[2]), int(pp[3]))

    print(picks)

    generate_test(
        args.test_dir, args.topology_file, 
        req_sizes, bf_sizes, entry_size_records, table_sizes, modes, picks)