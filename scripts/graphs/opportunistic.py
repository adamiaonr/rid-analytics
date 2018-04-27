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
from plot_utils import *

matplotlib.rcParams.update({'font.size': 16})
topology_labels = {1221: 'Telstra (1221)', 3356: 'Level3 (3356)', 3257: 'Tiscali (3257)', 4755: 'VSNL (4755)', 7018: 'AT&T (7018)'}

def add_label_cols(file_label, data):
    data['topology'] = int(file_label.split("-")[0])
    data['bf-size']  = int(file_label.split("-")[1])
    data['req-size'] = int(file_label.split("-")[2])
    data['fwd-size'] = int(file_label.split("-")[3])
    data['tab-size'] = int(file_label.split("-")[4])
    data['mode'] = int(file_label.split("-")[5])
    data['src:dst']  = str(file_label.split("-")[7] + ':' + file_label.split("-")[8])

def read_data(data_dir, file_types = ['outcomes'], filters = {}):

    # aggregate everything in a dataframe
    # FIXME: this will get big, fast...
    data = defaultdict(pd.DataFrame)
    for file_name in sorted(glob.glob(os.path.join(data_dir, '*.tsv'))):

        # extract file type
        # if not an outcome file, skip
        file_type = file_name.split("/")[-1].split(".")[0]
        if file_type not in file_types:
            continue

        # extract test label
        file_label = file_name.split(".")[1]
        # FIXME: comment this after testing
        if file_label.split("-")[0] not in ['1221', '3257', '3356', '7018']:
            continue 

        # apply filters
        if 'topology' in filters:
            if file_label.split("-")[0] != filters['topology']:
                continue

        if 'mode' in filters:
            if file_label.split("-")[5] != str(filters['mode']):
                continue

        # if 'case' in filters:
        #     if file_label.split("-")[9] != filters['case']:
        #         continue

        if 'bf-size' in filters:
            if file_label.split("-")[1] != filters['bf-size']:
                continue

        # extract data from file
        _data = pd.read_csv(file_name, sep = "\t").convert_objects(convert_numeric = True)
        # add label columns to _data
        add_label_cols(file_label, _data)

        # calculate the average latency per each outcome type
        #   - sum prob grouped by ['src:dst', 'type', 'latency']
        __data = _data.groupby(['src:dst', 'type', 'latency'])['prob'].agg('sum').reset_index()
        __data['avg-latency'] = __data['latency'] * __data['prob']
        #   - sum prob grouped by ['src:dst', 'type']
        _data = _data.groupby(['src:dst', 'type'])['prob'].agg('sum').reset_index()
        #   - add avg. latency to _data
        _data['avg-latency'] = __data.groupby(['src:dst', 'type'])['avg-latency'].agg('sum').reset_index()['avg-latency'] / _data['prob']

        # update main dataframe
        data[file_type] = pd.concat([data[file_type], _data], ignore_index = True)

    # calc sum of each diff. outcome (we use the term 'event' in the paper),
    # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size'} 
    for file_type in file_types:
        if file_type == 'outcomes':
            data[file_type] = data[file_type].groupby(['type'])[['prob', 'avg-latency']].agg('mean').reset_index()
        elif file_type == 'events':
            data[file_type] = data[file_type].groupby(['topology', 'bf-size', 'req-size', 'tab-size', 'fwd-size'])[['llm', 'mlm', 'slm']].agg('mean').reset_index()

    return data

def plot(data_dir, test_dir, output_dir):

    # labels & stuff
    topology_keys = ['1221', '3257', '3356', '7018']
    # grayscale for topologies
    topology_colors = ['#000000', '#708090', '#bebebe', 'white']
    # keys
    modes = [0, 3]
    cases = ['no-tps', 'tps']

    # prepare the graph parameters
    bar_width = 0.10

    # use the classic plot style
    plt.style.use('classic')

    # fig
    fig = plt.figure(figsize=(5, 3.5))
    # (1) avg. nr. of used links (bar chart)
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(False)
    # (2) fwd efficiency (line chart)
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, ls = 'dotted', lw = 0.50)

    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    pos = 0.0
    setlabel = True
    xbf = defaultdict()

    mstr = {'0192' : 'AML\n(192)', '0384' : 'AML\n(384 bit)', '3192' : 'FB\n(192 bit)'}
    cstr = {'no-tps' : 'No TPs', 'tps' : 'w/ TPs'}

    for c in cases:

        _data_dir = os.path.join(data_dir, ('%s/results' % (c)))

        for bf in ['384', '192']:
            for m in modes:

                if (m == 3) and (bf == '384'):
                    continue

                if (m == 0) and (bf == '192'):
                    continue

                for t, topology in enumerate(topology_keys):

                    _data = read_data(_data_dir, file_types = ['outcomes'], filters = {'topology' : topology, 'mode' : m, 'bf-size' : bf})['outcomes']
                    print("topology : %s, mode : %s, case : %s, bf-size : %s" % (t, m, c, bf))
                    print(_data)

                    label = ['', '', '']
                    if setlabel == True:
                        label[0] = topology_labels[int(topology)]
                        if t == 0:
                            label[1] = 'corr. del.'
                            label[2] = 'incorr. del.'

                    # bars hold avg. latency
                    print(_data.loc[_data['type'] == 'cd']['avg-latency'])
                    cd_lat = float(_data.loc[_data['type'] == 'cd']['avg-latency'])

                    ax2.bar(
                        pos, 
                        cd_lat,
                        color = topology_colors[t], linewidth = 0.5, width = bar_width, label = label[0])

                    # lines hold avg. # of deliveries
                    ax1.plot(
                        pos + (bar_width / 2.0),
                        _data.loc[_data['type'] == 'cd']['prob'],
                        linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

                    ax1.plot(
                        pos + (bar_width / 2.0),
                        _data.loc[_data['type'] == 'id']['prob'],
                        linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'v', label = label[2])

                    if t == 2:
                        xbf[pos] = str(mstr[str(m) + bf])

                    if (t == 3) and (m == 0):
                        xbf[pos + (1.25 * bar_width)] = '\n\n' + str(cstr[c])

                    pos += bar_width

                # only include lable once
                setlabel = False
                pos += (bar_width / 2.0)

        pos += (bar_width * 2.0)

    # an horizontal shaded area, for fwd efficiency
    ax1.axhspan(0, 2, linewidth = 0.0, facecolor = '#bebebe', alpha=0.50)

    # legend
    ax2.legend(fontsize = 10, ncol = 1, loc = 'upper right')
    ax1.legend(fontsize = 10, ncol = 1, loc = 'upper left')

    # axis labels
    ax1.set_xlabel("forwarding mode\nlocal annc. size")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax2.set_ylabel("avg. latency (# hops)")
    ax1.set_ylabel("avg. # deliveries")

    # avg used link axis
    ax2.set_xlim(0.0 - bar_width, pos - (1.5 * bar_width))
    ax2.set_ylim(0, 12)
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    
    # fwd effic axis
    ax1.set_ylim(-1, 5)
    ax1.set_yticks([0, 1, 2, 3])

    plt.savefig(os.path.join(output_dir, "opportunistic.pdf"), bbox_inches='tight', format = 'pdf')
