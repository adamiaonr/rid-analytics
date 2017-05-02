import argparse
import re
import sys
import glob
import os
import math
import binascii
import random
import subprocess
import time

import numpy as np
import networkx as nx
import matplotlib.pyplot as plt
import xml.etree.cElementTree as et

from xml.dom import minidom
from prettytable import PrettyTable

from collections import defaultdict
from collections import OrderedDict

def run_analysis(params):

    args = ("/home/adamiaonr/workbench/rid-analytics/run-analysis", 
        "--scn-file", params['scn-file'], 
        "--data-dir", params['data-dir'], 
        "--output-label", params['output-label'], 
        "--bf-size", params['bf-size'], 
        "--request-size", params['request-size'], 
        "--mm-mode", params['mm-mode'], 
        "--resolution-mode", params['resolution-mode'], 
        "--origin-server", params['origin-server'], 
        "--start-router", params['start-router'])

    start_time = time.time()
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(args, stdout=devnull)
    print("run_analysis for %s finished. [time : %s]\n" % (params['output-label'], time.time() - start_time))

    # command = "/home/adamiaonr/workbench/rid-analytics/run-analysis "
    # for arg in params:
    #     command += "--" + arg + " " + params[arg] + " "
    # command += "&> /dev/null"

    # start_time = time.time()
    # os.system(command)

    # popen = subprocess.Popen(args, stdout=subprocess.PIPE)
    # popen.wait()

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--test-file", 
         help = """path w/ topology files (e.g. .cch for isp maps, edge folders 
            for pop-level maps, etc.)""")

    args = parser.parse_args()

    if not args.test_file:
        sys.stderr.write("""%s: [ERROR] please supply a .test file path\n""" % sys.argv[0])
        parser.print_help()
        sys.exit(1)

    # extract info from .test file
    test_run = et.parse(args.test_file)
    test_run_root = test_run.getroot()

    params = defaultdict()
    for test in test_run.findall('test'):

        # test id to use as prefix to output labels
        test_id = test.get('id')

        params['data-dir'] = test.find('results_dir').text
        params['bf-size'] = test.find('bf-size').text
        params['request-size'] = test.find('req-size').text
        params['mm-mode'] = test.find('modes').text.split("-", 1)[0]
        params['resolution-mode'] = test.find('modes').text.split("-", 1)[1]

        for path in test.find('paths').findall('path'):

            params['scn-file'] = path.get('file')
            params['start-router'] = str(int(path.text.split(",")[0]))
            params['origin-server'] = str(int(path.text.split(",")[-1]))

            length = path.get('length')
            outdegree = path.get('outdegre')
            params['output-label'] = ("%s-%s-%s-%s-%s" % ( 
                test_id, length, outdegree,
                params['start-router'], params['origin-server']))

            run_analysis(params)
