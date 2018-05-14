import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
import sys
import glob

from datetime import date
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict

# custom imports
import plot_utils
import globalr
import cdn
import opportunistic
import tradeoff
import event_breakdown

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--data-dir", 
         help = """dir(s) w/ .tsv files.""")
    parser.add_argument(
        "--output-dir", 
         help = """dir on which to print graphs.""")
    parser.add_argument(
        "--case", 
         help = """the case you want to output. e.g. 'base'.""")
    parser.add_argument(
        "--subcase", 
         help = """the sub-case you want to output. e.g. 'bfs'.""")
    parser.add_argument(
        "--test-file", 
         help = """.test file path""")

    args = parser.parse_args()

    # quit if a dir w/ causality files hasn't been provided
    if not args.data_dir:
        sys.stderr.write("""%s: [ERROR] please supply a data dir!\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    if not args.output_dir:
        sys.stderr.write("""%s: [ERROR] please supply an output dir!\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    if not args.test_file:
        sys.stderr.write("""%s: [ERROR] please supply a .test file path!\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    if args.case == 'global':

        if not args.subcase:
            sys.stderr.write("""%s: [ERROR] please supply a subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

        if args.subcase == 'bf-sizes':
            globalr.plot_bf_sizes(args.data_dir, args.test_file, args.output_dir)
        elif args.subcase == 'tab-sizes':
            globalr.plot_tab_sizes(args.data_dir, args.test_file, args.output_dir)
        elif args.subcase == 'req-sizes':
            globalr.plot_req_sizes(args.data_dir, args.test_file, args.output_dir)
        elif args.subcase == 'tab-req-size-tradeoff':
            globalr.plot_table_req_size_tradeoff(args.data_dir, args.test_file, args.output_dir)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

    elif args.case == 'cdn':

        if not args.subcase:
            sys.stderr.write("""%s: [ERROR] please supply a subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

        if args.subcase == 'usage':
            cdn.plot_usage(args.data_dir, args.test_file, args.case, args.output_dir)
        elif args.subcase == 'breakdown':
            cdn.plot_breakdown(args.data_dir, args.test_file, args.case, args.output_dir)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

    elif args.case == 'opportunistic':

        if args.subcase == 'usage':
            cdn.plot_usage_(args.data_dir, args.test_file, args.case, args.output_dir)
        elif args.subcase == 'breakdown':
            cdn.plot_breakdown_(args.data_dir, args.test_file, args.case, args.output_dir)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

    elif args.case == 'cdn-lpm':

        if args.subcase == 'usage':
            cdn.plot_usage_(args.data_dir, args.test_file, args.case, args.output_dir)
        elif args.subcase == 'breakdown':
            cdn.plot_breakdown_(args.data_dir, args.test_file, args.case, args.output_dir)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

    elif args.case == 'tradeoff':
        tradeoff.plot_tradeoff(args.data_dir, args.test_file, args.output_dir)

    else:
        sys.stderr.write("""%s: [ERROR] please supply a valid case\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)