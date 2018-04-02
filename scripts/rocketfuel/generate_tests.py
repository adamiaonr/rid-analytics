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

def add_new_test(
    main_block, 
    topology_obj, 
    test_parameters, 
    path_examples, 
    test_dir):

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
    for path_key in path_examples:

        # save the scenario info about the particular test
        path = path_examples[path_key]
        
        # generate the .scn filename for this scenario
        scn_filename = ""
        if len(add_suffixes) > 0:
            scn_filename_suffix = '-'.join([x for x in add_suffixes])
            scn_filename = ("configs/%s-%03d-%03d-%s.scn" % (test_id, path[0], path[-1], scn_filename_suffix))
        else:
            scn_filename = ("configs/%s-%03d-%03d.scn" % (test_id, path[0], path[-1]))                
        scn_filename = os.path.join(test_dir, scn_filename)

        new_path_block = et.SubElement(paths_block, "path", length = str(len(path) - 1), avg_outdegree = path_key.split(":")[-1], file = scn_filename).text = str(','.join(("%03d" % (x)) for x in path))

        # add true positives
        _topology = topology_obj.topology.copy()
        _topology_obj = Topology(_topology)
        _topology_obj.set_shortest_paths(topology_obj.get_shortest_paths())

        _topology_obj.add_tp_route(path[-1], int([x for x in test_parameters['entry-sizes']][0]))
        print("%s::add_new_test() : [INFO] added tp route : %s" % (sys.argv[0], path))

        # add a secondary tp source, if required (??)
        # if len(tps) > 0:
        #     for tp in tps:
        #         _topology_obj.add_tp_route(
        #             tp_src_id = -1, 
        #             tp_size = int(tp.split(":")[1]), 
        #             radius = int(tp.split(":")[2]), 
        #             n_tp_srcs = int(tp.split(":")[0]), 
        #             dst_id = path[0])

        # add fp records
        for fp_record in test_parameters['fps']:

            fp_record = fp_record.split(":")
            if fp_record[0] == 'S':

                _topology_obj.add_fp_route(
                    src_id = path[0],
                    fp_size = int(fp_record[1]),
                    fp_size_proportion = int(fp_record[2]),
                    radius = int(fp_record[3]))

            else:
                print("not adding additional fp sources")

        # now save the actual .scn file
        _topology_obj.convert_to_scn(
            test_parameters['entry-sizes'], 
            test_parameters['table-size'], 
            test_parameters['req-size'], 
            scn_filename)

    results_dir_block = et.SubElement(test_block, "results_dir").text = os.path.join(test_dir, "results/")

    return 0

def generate_test(
    test_dir, 
    topology_file, 
    req_sizes, 
    bf_sizes, 
    entry_size_records, 
    table_sizes, 
    fps, 
    tps, 
    modes, 
    path_sizes, 
    selected_paths,
    add_suffixes):

    # describe the test to be generated
    print("%s::generate_test() : [INFO] test info:" % (sys.argv[0]))
    print("\ttest dir: %s" % (str(test_dir)))
    print("\ttopology file: %s" % (str(topology_file)))
    print("\treq sizes: %s" % (str(req_sizes)))
    print("\tbf sizes: %s" % (str(bf_sizes)))
    print("\tentry sizes: %s" % (str(entry_size_records)))
    print("\ttable sizes: %s" % (str(table_sizes)))
    print("\tfalse pos. entries: %s" % (str(fps)))
    print("\ttrue pos. entries: %s" % (str(tps)))
    print("\tmodes: %s" % (str(modes)))
    print("\tpath sizes: %s" % (str(path_sizes)))
    print("\tselected paths: %s" % (str(selected_paths)))
    print("\tsuffixes: %s" % (str(add_suffixes)))

    # if test_dir doesn't contain 3 main folders, create them now
    for folder in ['configs', 'results', 'topologies']:
        p = os.path.join(test_dir, folder)
        if not os.path.exists(p):
            os.makedirs(p)

    main_block = et.Element("test_run")

    # extract the topology nr.
    topology_nr = int(topology_file.split("/")[-2])
    topology_nr_block = et.SubElement(main_block, "topology", file = topology_file).text = ("%d" % (topology_nr))

    # build the basic topology (no .scn file yet), and generate path examples
    topology_obj = parse_pop_level_map(topology_file)
    topology_obj.draw_pop_level_map(os.path.join(test_dir, ("topologies/%d.pdf" % (topology_nr))))

    # get n examples of paths of size s
    if len(selected_paths) > 0:
        path_examples = topology_obj.get_paths(selected_paths)
    else:
        path_examples = topology_obj.generate_paths(int(path_sizes.split(":")[0]), int(path_sizes.split(":")[1]))

    print("%s::generate_test() : [INFO] generated %d path examples (n : %d, size : %d, topology : %d):" 
        % (sys.argv[0], len(path_examples), int(path_sizes.split(":")[0]), int(path_sizes.split(":")[1]), topology_nr))
    print(path_examples)

    # build the test cases w/ good ol' nested for loops
    for req_size in req_sizes:
        for bf_size in bf_sizes:
            for record in entry_size_records:

                entry_sizes = defaultdict();
                entry_size_proportions = record.split("|")

                for p in entry_size_proportions:
                    entry_sizes[p.split(":", 1)[0]] = (float(p.split(":", 1)[1]) / 100.0)

                for table_size in table_sizes:
                    for mode in modes:

                        test_parameters = defaultdict()

                        test_parameters['topology-nr'] = topology_nr
                        test_parameters['req-size'] = req_size
                        test_parameters['bf-size'] = bf_size
                        test_parameters['entry-sizes'] = entry_sizes
                        test_parameters['table-size'] = int(table_size)
                        test_parameters['fps'] = fps
                        test_parameters['tps'] = tps
                        test_parameters['modes'] = mode.replace(":", "-")
                        test_parameters['add-suffixes'] = add_suffixes

                        add_new_test(main_block, topology_obj, test_parameters, path_examples, test_dir)

    xmlstr = minidom.parseString(et.tostring(main_block)).toprettyxml(indent="    ")
    with open(os.path.join(test_dir, ("configs/%d.test" % (topology_nr))), "w") as f:
        f.write(xmlstr)

    return 0

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--topology-file", 
         help = """filename of edge.wt topology file (or dir for 
            'print-stats' option)""")

    parser.add_argument(
        "--print-stats", 
         help = """print statistics about RocketFuel dataset""",
         action = "store_true")

    # options (self-explanatory)
    parser.add_argument(
        "--test-dir", 
         help = """dir on which to save .test and .scn files""")

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
        "--add-tps", 
         help = """add tp sources. syntax is 
            <num_tp_srcs>:<annc. size>:<annc. radius>. e.g. '--add-tps "1:2:2' """)

    parser.add_argument(
        "--add-fps", 
         help = """add fp announcements around some node n. syntax 
            is <n id>:<annc. size>:<size %%>:<annc. radius>. if <n id> == 'S' 
            sets <n id> to the source of request. e.g. '--add-fp "S:2:50:2' """)

    parser.add_argument(
        "--modes", 
         help = """MM_MODE (0 for 'flood', 1 for 'random', or 2 for 'fallback'), 
            RES_MODE (0 for 'drop packets', 1 for 'resolve w/ fallback'), in that order, 
            separated by ':' and '|' e.g. '--modes 0:0|0:1'""")

    parser.add_argument(
        "--path-sizes", 
         help = """pick (up to) n paths of a particular size s. 
            e.g. '--path-sizes 10:4' for (up to) 10 paths of size 4.""")

    parser.add_argument(
        "--selected-paths", 
         help = """specify the paths to choose for scenarios (instead of randomly 
            generating them). 
            e.g. '--selected-paths 0:3|0:2' for 2 paths starting at node 0, and 
            ending at nodes 3 and 2, respectively.""")

    parser.add_argument(
        "--add-suffix", 
         help = """add an extra suffix to the results and .scn file.""")

    args = parser.parse_args()

    if args.print_stats:
        print_pop_level_statistics(args.topology_file)
        sys.exit(0)

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

    fps = []
    if args.add_fps:
        fps = args.add_fps.split("|")

    tps = []
    if args.add_tps:
        tps = args.add_tps.split("|")

    modes = []
    if args.modes:
        modes = args.modes.split("|")

    if not args.path_sizes:
        args.path_sizes = "1:4"

    selected_paths = []
    if args.selected_paths:
        selected_paths = args.selected_paths.split("|")

    add_suffixes = []
    if args.add_suffix:
        add_suffixes = args.add_suffix.split(":")

    generate_test(
        args.test_dir, args.topology_file, 
        req_sizes, bf_sizes, entry_size_records, table_sizes, fps, tps, modes, args.path_sizes, selected_paths,
        add_suffixes)