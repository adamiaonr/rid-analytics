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
import multiprocessing as mp

from xml.dom import minidom
from prettytable import PrettyTable

from collections import defaultdict
from collections import OrderedDict

def run_analysis(
    scn_file,
    data_dir,
    output_label,
    bf_size,
    request_size,
    mm_mode,
    resolution_mode,
    origin_server,
    start_router):

    args = ("/home/adamiaonr/workbench/rid-analytics/run-analysis", 
        "--scn-file", scn_file, 
        "--data-dir", data_dir, 
        "--output-label", output_label, 
        "--bf-size", bf_size, 
        "--request-size", request_size, 
        "--mm-mode", mm_mode, 
        "--resolution-mode", resolution_mode, 
        "--origin-server", origin_server, 
        "--start-router", start_router)

    print(args)

    start_time = time.time()
    with open(os.devnull, 'w') as devnull:
        subprocess.check_call(args, stdout=devnull)
    
    print("run_analysis for %s finished. [time : %s]\n" % (output_label, time.time() - start_time))

    return output_label

if __name__ == "__main__":

    # use an ArgumentParser for a nice CLI
    parser = argparse.ArgumentParser()

    # options (self-explanatory)
    parser.add_argument(
        "--test-file", 
         help = """.test file w/ description of tests to run""")

    args = parser.parse_args()

    if not args.test_file:
        sys.stderr.write("""%s: [ERROR] please supply a .test file path\n""" % sys.argv[0])
        parser.print_help()
        sys.exit(1)

    # extract info from .test file
    test_run = et.parse(args.test_file)
    test_run_root = test_run.getroot()

    params = dict()
    for test in test_run.findall('test'):

        # test id to use as prefix to output labels
        test_id = test.get('id')

        params['data-dir'] = test.find('results_dir').text
        params['bf-size'] = test.find('bf-size').text
        params['request-size'] = test.find('req-size').text
        params['mm-mode'] = test.find('modes').text.split("-", 1)[0]
        params['resolution-mode'] = test.find('modes').text.split("-", 1)[1]

        # we use a thread pool to run tests in parallel
        pool = mp.Pool(mp.cpu_count())
        tasks = []

        for path in test.find('paths').findall('path'):

            params['scn-file'] = path.get('file')
            params['start-router'] = str(int(path.text.split(",")[0]))
            params['origin-server'] = str(int(path.text.split(",")[-1]))

            params['output-label'] = ("%s-%s-%s" % ( 
                test_id,
                params['start-router'], params['origin-server']))

            tasks.append((
                params['scn-file'], 
                params['data-dir'], 
                params['output-label'], 
                params['bf-size'], 
                params['request-size'], 
                params['mm-mode'], 
                params['resolution-mode'], 
                params['origin-server'], 
                params['start-router']))

        if len(tasks) > 0:

            jobs_remaining = len(tasks)
            results = [pool.apply_async(run_analysis, task) for task in tasks]

            for result in results:

                jobs_remaining = jobs_remaining - 1
                test_label = result.get()
                
                if test_label is not None:
                    print("finished %s. %d jobs remaining." 
                        % (test_label, jobs_remaining))

        # keep things tidy
        pool.close()
        pool.join()