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
topology_labels = {1221: 'Telstra', 3356: 'Level3', 3257: 'Tiscali', 4755: 'VSNL', 7018: 'AT&T', 1337 : 'Ideal'}

def add_label_cols(file_label, data):
    data['topology'] = int(file_label.split("-")[0])
    data['bf-size']  = int(file_label.split("-")[1])
    data['req-size'] = int(file_label.split("-")[2])
    data['fwd-size'] = int(file_label.split("-")[3])
    data['tab-size'] = int(file_label.split("-")[4])
    data['mode'] = int(file_label.split("-")[5])
    data['src:dst']  = str(file_label.split("-")[7] + ':' + file_label.split("-")[8])
    data['case'] = str(file_label.split("-")[9])

def read_data(data_dir, file_types = ['outcomes'], filters = {}):

    print(file_types)

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

        if 'case' in filters:
            if file_label.split("-")[9] != filters['case']:
                continue

        if 'bf-size' in filters:
            if file_label.split("-")[1] != filters['bf-size']:
                continue

        # extract data from file
        _data = pd.read_csv(file_name, sep = "\t").convert_objects(convert_numeric = True)
        # add label columns to _data
        add_label_cols(file_label, _data)

        if file_type == 'outcomes':
            # calculate the average latency per each outcome type
            #   - sum prob grouped by ['src:dst', 'type', 'latency']
            __data = _data.groupby(['src:dst', 'type', 'latency'])['prob'].agg('sum').reset_index()
            __data['avg-latency'] = __data['latency'] * __data['prob']
            #   - sum prob grouped by ['src:dst', 'type']
            _data = _data.groupby(['src:dst', 'type'])['prob'].agg('sum').reset_index()
            #   - add avg. latency to _data
            _data['avg-latency'] = __data.groupby(['src:dst', 'type'])['avg-latency'].agg('sum').reset_index()['avg-latency'] / _data['prob']

        elif file_type == 'events':
            # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size', 'src:dst'}, and sum() the diff event types
            _data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst'])['llm', 'mlm', 'slm'].agg('sum').reset_index()

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

def plot_no_tps(data_dir, test_dir, output_dir):

    # labels & stuff
    # topology_keys = ['1221', '3257', '3356', '7018']
    topology_keys = ['1221', '3257', '3356', '7018', '1337']
    # grayscale for topologies
    topology_colors = ['#000000', '#708090', '#bebebe', 'white']
    # keys
    modes = [0, 3]
    bf_sizes = {0 : 384, 3 : 192}
    cases = ['S212', 'S512']

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

    mstr = {'0192' : 'AML\n(192)', '0384' : 'AML (384 bit)', '3192' : 'FB (192 bit)'}
    cstr = {'S212' : 2, 'S512' : 5}

    for t, topology in enumerate(topology_keys):
        for c in cases:
            for m in modes:

                cd_lat = 0.0
                link_usage = 0.0

                if topology != '1337':

                    _data = read_data(data_dir, file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : c, 'mode' : m})
                    # bars  : link usage
                    # lines : avg. latency
                    cd_lat = float(_data['outcomes'].loc[_data['outcomes']['type'] == 'cd']['avg-latency'])

                    link_usage = _data['events']['mlm'].values + _data['events']['slm'].values
                    if m == 3:
                        __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']
                        link_usage = link_usage + __d['prob'].values * (4.0 - __d['avg-latency'].values)

                else:
                    cd_lat = 4.0
                    link_usage = 4.0

                bf = bf_sizes[m]
                print("topology : %s, mode : %s, case : %s, bf-size : %s" % (t, m, c, bf))

                label = ['', '', '']
                if setlabel == True:
                    label[0] = str(mstr[str(m) + str(bf)])
                    if m == 0:
                        label[1] = 'avg. latency'

                ax2.bar(
                    pos, 
                    link_usage,
                    color = topology_colors[m], linewidth = 0.5, width = bar_width, label = label[0])

                # lines hold avg. # of deliveries
                ax1.plot(
                    pos + (bar_width / 2.0),
                    cd_lat,
                    linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

                # if t == 2:
                #     xbf[pos] = str(mstr[str(m) + bf])

                # if (t == 3) and (m == 0):
                #     xbf[pos + (1.25 * bar_width)] = '\n\n' + str(cstr[c])

                if m == 0:
                    xbf[pos + (bar_width)] = str(cstr[c])

                pos += bar_width

            if c == 'S212':
                xbf[pos + (0.25 * bar_width)] = '\n' + str(topology_labels[int(topology)])

            # only include lable once
            setlabel = False
            pos += (bar_width / 2.0)

        pos += (bar_width * 1.0)

    # an horizontal shaded area, for fwd efficiency
    ax1.axhspan(3, 5, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.25)

    ax1.plot(
        [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        [3.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    ax1.plot(
        [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        [5.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    # legends
    leg = []
    leg.append(ax1.legend(fontsize = 10, ncol = 1, loc = 'upper left'))
    leg.append(ax2.legend(fontsize = 10, ncol = 1, loc = 'upper right', title = 'fwd. strategy'))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 10)

    # axis labels
    ax1.set_xlabel("cache annc. size\ntopology")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax2.set_ylabel("avg. # of used links")
    ax1.set_ylabel("avg. latency")

    # avg used link axis
    ax2.set_xlim(0.0 - (bar_width / 2.0), pos - (1.0 * bar_width))
    ax2.set_ylim(0, 10)
    ax2.set_yticks([0, 1, 2, 3, 4, 5, 6])
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    # ax2.set_ylim(1, 10000)
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 5, 1)])
    # ax2.tick_params(axis = 'y', which = 'minor', bottom = 'off')
    
    # fwd effic axis
    ax1.set_ylim(-2, 8)
    ax1.set_yticks([3.0, 4.0, 5.0])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "cdn-no-tps.pdf"), format = 'pdf')

def plot_tps(data_dir, test_dir, output_dir):

    # labels & stuff
    topology_keys = ['1221', '3257', '3356', '7018']
    # grayscale for topologies
    topology_colors = ['#000000', '#708090', '#bebebe', 'white']
    # keys
    modes = [0, 3]
    cases = ['S212', 'S512']

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
    cstr = {'S212' : 2, 'S512' : 5}

    for c in cases:
        for bf in ['384', '192']:
            for m in modes:

                if (m == 3) and (bf == '384'):
                    continue

                if (m == 0) and (bf == '192'):
                    continue

                for t, topology in enumerate(topology_keys):

                    _data = read_data(data_dir, file_types = ['outcomes'], filters = {'topology' : topology, 'case' : c, 'mode' : m, 'bf-size' : bf})['outcomes']
                    print("topology : %s, mode : %s, case : %s, bf-size : %s" % (t, m, c, bf))
                    print(_data)

                    label = ['', '', '']
                    if setlabel == True:
                        label[0] = topology_labels[int(topology)]
                        if t == 0:
                            label[1] = 'corr. del.'
                            label[2] = 'incorr. del.'

                    # bars hold avg. latency
                    cd_lat      = float(_data.loc[_data['type'] == 'cd']['avg-latency'])

                    ax2.bar(
                        pos, 
                        cd_lat,
                        color = topology_colors[t], linewidth = 0.5, width = bar_width, label = label[0])

                    # lines hold avg. # of deliveries
                    ax1.plot(
                        pos + (bar_width / 2.0),
                        _data.loc[_data['type'] == 'cd']['prob'],
                        linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

                    id_del = 0.0
                    df = _data.loc[_data['type'] == 'id']
                    if not df.empty:
                        id_del = float(df['prob'])
                    ax1.plot(
                        pos + (bar_width / 2.0),
                        id_del,
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
    ax2.set_ylim(0, 10)
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    
    # fwd effic axis
    ax1.set_ylim(-1, 4)
    ax1.set_yticks([0, 1, 2])

    plt.savefig(os.path.join(output_dir, "cdn-tps.pdf"), bbox_inches='tight', format = 'pdf')