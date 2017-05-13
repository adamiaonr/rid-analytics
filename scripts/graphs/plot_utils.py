import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
import sys
import glob

import xml.etree.cElementTree as et

from datetime import date
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict
from itertools import chain, izip

# interface event codes
EVENT_NLM = 0 # no link matches
EVENT_MLM = 1 # multiple link matches
EVENT_LLM = 2 # local link match
EVENT_SLM = 3 # single link match (other than local)
EVENT_TTL = 4 # drop due to rtt expiration

# outcome codes
OUTCOME_CORRECT_DELIVERY     = 0
OUTCOME_INCORRECT_DELIVERY   = 1
OUTCOME_FALLBACK_DELIVERY    = 2
OUTCOME_FALLBACK_RELAY       = 3
OUTCOME_PACKET_DROP          = 4
OUTCOME_TTL_DROP             = 5
OUTCOME_UNDEF                = 6

ileave = lambda *iters: list(chain(*izip(*iters)))

# full method with doctests
def interleave_n(*iters):
    """
    Given two or more iterables, return a list containing 
    the elements of the input list interleaved.
    
    >>> x = [1, 2, 3, 4]
    >>> y = ('a', 'b', 'c', 'd')
    >>> interleave(x, x)
    [1, 1, 2, 2, 3, 3, 4, 4]
    >>> interleave(x, y, x)
    [1, 'a', 1, 2, 'b', 2, 3, 'c', 3, 4, 'd', 4]
    
    On a list of lists:
    >>> interleave(*[x, x])
    [1, 1, 2, 2, 3, 3, 4, 4]
    
    Note that inputs of different lengths will cause the 
    result to be truncated at the length of the shortest iterable.
    >>> z = [9, 8, 7]
    >>> interleave(x, z)
    [1, 9, 2, 8, 3, 7]
    
    On single iterable, or nothing:
    >>> interleave(x)
    [1, 2, 3, 4]
    >>> interleave()
    []
    """
    return list(chain(*izip(*iters)))

def interleave(a, b):
    c = list(zip(a, b))
    return [elt for sublist in c for elt in sublist]

def extract_data(data_dir):

    data = defaultdict(OrderedDict)

    for file_name in sorted(glob.glob(os.path.join(data_dir, '*.tsv'))):

        file_type = file_name.split(".")[0].split("/")[-1]
        file_label = file_name.split(".")[1]

        # print("filename = %s, type = %s, label = %s" % (file_name, file_type, file_label))

        data[file_type][file_label] = pd.read_csv(file_name, sep = "\t")
        data[file_type][file_label] = data[file_type][file_label].convert_objects(convert_numeric = True)

    return data

def get_path(test_file, file_label):

    print(test_file)
    print(file_label)

    # extract info from .test file
    test_run = et.parse(test_file)
    test_run_root = test_run.getroot()

    path_lengths = defaultdict()
    avg_outdegrees = defaultdict()

    # use test_id to get additional information about the 
    # test from the .test file
    test_id_components = file_label.split("-")

    src_id = 0
    dst_id = 0
    if ('tp' in file_label) or ('fp' in file_label):
        src_id = int(test_id_components[8])
        dst_id = int(test_id_components[9])
    else:
        src_id = int(test_id_components[7])
        dst_id = int(test_id_components[8])

    test_id = ''
    if 'tp5' in file_label:
        test_id = '-'.join(test_id_components[:7])
        test_id += '-tp5'
        print(test_id)
    elif 'fp5' in file_label:
        test_id = '-'.join(test_id_components[:7])
        test_id += '-fp5'
    elif 'tp' in file_label:
        test_id = '-'.join(test_id_components[:8])
    else:
        test_id = '-'.join(test_id_components[:7])

    # print("get_path_length() : test_id = %s" % (test_id))

    print(test_id)

    for test in test_run.findall('test'):
        # test id to use as prefix to output labels
        if test.get('id') != test_id:
            continue

        for path in test.find('paths').findall('path'):
            if int(path.text.split(',')[0]) == src_id and int(path.text.split(',')[-1]) == dst_id:
                return [int(p) for p in path.text.split(',')]

def get_path_info(test_file, file_label):

    # extract info from .test file
    test_run = et.parse(test_file)
    test_run_root = test_run.getroot()

    path_lengths = defaultdict()
    avg_outdegrees = defaultdict()

    # use test_id to get additional information about the 
    # test from the .test file
    test_id_components = file_label.split("-")
    if 'R' in file_label:
        del test_id_components[7]
    test_id = '-'.join(test_id_components[:7])
    # print("get_path_length() : test_id = %s" % (test_id))

    for test in test_run.findall('test'):
        # test id to use as prefix to output labels
        if test.get('id') != test_id:
            continue

        for path in test.find('paths').findall('path'):

            length_key = path.get('length')
            outdegree_key = path.get('outdegre')

            if length_key not in path_lengths:
                path_lengths[length_key] = defaultdict()
                avg_outdegrees[length_key] = defaultdict()

            path_lengths[length_key][outdegree_key] = (len(path.text.split(',')) - 1)
            avg_outdegrees[length_key][outdegree_key] = float(path.get('avg_outdegree'))

    return path_lengths, avg_outdegrees
