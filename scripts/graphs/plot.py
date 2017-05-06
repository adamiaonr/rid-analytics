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
import plot_base
import plot_global
import plot_utils

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--data-dir", 
         help = """dir w/ .tsv files.""")
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

    # if an output dir is not specified, use data-dir
    if not args.output_dir:
        args.output_dir = args.data_dir

    if not args.test_file:
        sys.stderr.write("""%s: [ERROR] please supply a .test file path!\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    # extract the data from all files in the data dir
    data = plot_utils.extract_data(args.data_dir)

    if args.case == 'base':
        plot_base.plot_efficiency(data, args.test_file)

    elif args.case == 'global':

        if not args.subcase:
            sys.stderr.write("""%s: [ERROR] please supply a subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

        if args.subcase == 'bf-sizes':
            plot_global.plot_bf_sizes(data, args.test_file)
        elif args.subcase == 'request-sizes':
            plot_global.plot_request_sizes(data, args.test_file)
        elif args.subcase == 'bf-sizes-fallbacks':
            plot_global.plot_bf_sizes_fallbacks(data, args.test_file)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase for %s case\n""" % (sys.argv[0], args.case))
            parser.print_help()
            sys.exit(1)

    else:
        sys.stderr.write("""%s: [ERROR] please supply a valid case\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)