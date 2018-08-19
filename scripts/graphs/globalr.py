import pandas as pd
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import os
import argparse
import sys
import glob
import math

from datetime import date
from datetime import datetime
from collections import defaultdict
from collections import OrderedDict

# custom imports
from plot_utils import *

matplotlib.rcParams.update({'font.size': 16})
# topology_labels = {1221: 'Telstra (1221)', 3356: 'Level3 (3356)', 3257: 'Tiscali (3257)', 4755: 'VSNL (4755)', 7018: 'AT&T (7018)'}
topology_labels = {1221: 'Telstra', 3356: 'Level3', 3257: 'Tiscali', 4755: 'VSNL', 7018: 'AT&T', 'ideal' : 'Ideal'}

def add_label_cols(file_label, data):
    # add columns w/ topology, bf size, req size and fwd size
    data['topology'] = int(file_label.split("-")[0])
    data['bf-size']  = int(file_label.split("-")[1])
    data['req-size'] = int(file_label.split("-")[2])
    data['fwd-size'] = int(file_label.split("-")[3])
    data['tab-size'] = int(file_label.split("-")[4])
    data['src:dst']  = str(file_label.split("-")[7] + ':' + file_label.split("-")[8])

def read_data(data_dir, file_types = ['outcomes', 'events']):

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

        # extract data from file
        _data = pd.read_csv(file_name, sep = "\t").convert_objects(convert_numeric = True)
        # add label columns to _data
        add_label_cols(file_label, _data)

        if file_type == 'outcomes':
            # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size', 'src:dst', 'type'}, and sum() by 'prob'
            _data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst', 'type'])['prob'].agg('sum').to_frame().reset_index()

        elif file_type == 'events':
            # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size', 'src:dst'}, and sum() the diff event types
            _data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst'])['llm', 'mlm', 'slmc', 'slmi', 'fpd'].agg('sum').reset_index()

        # update main dataframe
        data[file_type] = pd.concat([data[file_type], _data], ignore_index = True)

    # calc sum of each diff. outcome (we use the term 'event' in the paper),
    # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size'} 
    for file_type in file_types:
        if file_type == 'outcomes':
            data[file_type] = data[file_type].groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'type'])['prob'].agg('mean').reset_index()
        elif file_type == 'events':
            data[file_type] = data[file_type].groupby(['topology', 'bf-size', 'req-size', 'tab-size', 'fwd-size'])[['llm', 'mlm', 'slmc', 'slmi', 'fpd']].agg('mean').reset_index()

    return data

def plot_table_req_size_tradeoff(data_dir, test_dir, output_dir):
    
    data = read_data(data_dir, ['events'])

    tp_keys = sorted(list(data['events']['topology'].unique()))
    rs_keys = sorted(list(data['events']['req-size'].unique()))
    ts_keys = sorted(list(data['events']['tab-size'].unique()))

    linestyles = ['-','--']
    markers = ['o', '^', 'v', 's']

    # use the classic plot style
    plt.style.use('classic')
    # avoid Type 3 fonts
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # fig
    fig = plt.figure(figsize = (5, 3.5))
    xlims = [0.5, 0.5]
    ylims = [2000, 200]

    ax = [None, None]
    for tt, tp in enumerate(tp_keys):

        ax[tt] = fig.add_subplot(211 + tt)
        ax[tt].xaxis.grid(True)
        ax[tt].yaxis.grid(True)

        ax[tt].set_title("%s (%s)" % (topology_labels[tp], tp), fontsize = 12)

        pos = 0.0
        for t, ts in enumerate(ts_keys):

            _data = data['events'].loc[(data['events']['tab-size'] == ts) & (data['events']['topology'] == tp)]
            ts_series = []
            for rs in rs_keys:

                __data = _data.loc[_data['req-size'] == rs]
                ts_series.append((__data['mlm'].values + __data['slmi'].values + __data['slmc'].values) / 4.0)

            print(int(math.log10(int(ts))))
            print(np.arange(0, len(ts_series), 1))
            print(ts_series)
            ax[tt].plot(
                np.arange(3, len(ts_series) + 3, 1),
                ts_series,
                linewidth = 1.0, color = 'black', linestyle = '-', markersize = 5, marker = markers[t], 
                label = '10e' + str(int(math.log10(int(ts)))))

        ax[tt].set_ylabel("link usage\nexcess")
        ax[tt].set_yscale("log", nonposy = 'clip')
        ax[tt].set_ylim(xlims[tt], ylims[tt])

        ax[tt].set_xlim(-0.25 + 3, len(rs_keys) - 0.75 + 3)

    # ax[tt].plot(
    #     np.arange(5, len(ts_series) + 5, 1),
    #     [4.0] * len(ts_series),
    #     linewidth = 1.0, color = 'black', linestyle = '--')

    # legend
    leg = []
    leg.append(ax[0].legend(
        fontsize = 12, 
        ncol = 2, loc = 'upper left', title = 'table size', 
        handletextpad = 0.2, handlelength = 0.575, labelspacing = 0.2, columnspacing = 0.3,
        bbox_to_anchor=(0.0, 1.100)))
    leg.append(ax[1].legend(
        fontsize = 12, 
        ncol = 2, loc = 'upper left', title = 'table size', 
        handletextpad = 0.2, handlelength = 0.575, labelspacing = 0.2, columnspacing = 0.3,
        bbox_to_anchor=(0.0, 1.100)))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 12)

    # axis labels
    ax[1].set_xlabel("request size")
    # ax1.set_yticks([10**(x) for x in np.arange(0, 3, 1)])

    fig.subplots_adjust(left = None, bottom = None, right = None, top = None, wspace = None, hspace = 0.00)
    plt.tight_layout(h_pad = 0.10)
    plt.savefig(os.path.join(output_dir, "global-table-req-size-tradeoff.pdf"), format = 'pdf')


def plot_req_sizes(data_dir, test_dir, output_dir):

    # get data
    data = read_data(data_dir, ['events'])
    for file_type in ['events']:
        print(data[file_type])

    # labels & stuff
    topology_keys = sorted(list(data['events']['topology'].unique()))
    # grayscale for topologies
    topology_colors = ['#000000', '#708090', '#bebebe', 'white']
    rs_keys = sorted(list(data['events']['req-size'].unique()))

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
    for rs in rs_keys:
        for t, topology in enumerate(topology_keys):

            _data = data['events'].loc[
                (data['events']['topology'] == topology)
                & (data['events']['req-size'] == rs)]

            label = ['', '']
            if setlabel == True:
                label[0] = topology_labels[topology]
                if t == 0:
                    label[1] = 'fwd. effic.'

            ax2.bar(
                pos, 
                _data['llm'] + _data['slm'] - _data['mlm'],
                color = topology_colors[t], linewidth = 0.5, width = bar_width, label = label[0])

            ax1.plot(
                pos + (bar_width / 2.0),
                ((4.0 / (_data['llm'] + _data['slm'] - _data['mlm'])) * 2.0),
                linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

            if t == 2:
                xbf[pos] = str(rs)

            pos += bar_width

        # only include lable once
        setlabel = False
        pos += (bar_width * 2.0)

    # an horizontal shaded area, for fwd efficiency
    ax1.axhspan(0, 2, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    # legend
    ax2.legend(fontsize = 10, ncol = 1, loc = 'upper right')
    ax1.legend(fontsize = 10, ncol = 1, loc = 'upper left')

    # axis labels
    ax1.set_xlabel("request size\n\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax2.set_ylabel("avg. # of used links")
    ax1.set_ylabel("fwd. efficiency [0,1]")

    # avg used link axis
    ax2.set_xlim(0.0 - bar_width, pos - (1.0 * bar_width))
    ax2.set_ylim(1, 10000000)
    ax2.set_yscale("log", nonposy = 'clip')
    ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    
    # fwd effic axix
    ax1.set_ylim(-2, 5)
    ax1.set_yticks([0, 1, 2])
    ax1.set_yticklabels([0, 0.5, 1])

    plt.savefig(os.path.join(output_dir, "global-req-sizes.pdf"), bbox_inches='tight', format = 'pdf')

def plot_tab_sizes(data_dir, test_dir, output_dir):

    # get data
    data = read_data(data_dir, ['events'])
    for file_type in ['events']:
        print(data[file_type])

    # labels & stuff
    topology_keys = sorted(list(data['events']['topology'].unique()))
    # grayscale for topologies
    topology_colors = ['#000000', '#708090', '#bebebe', 'white']
    ts_keys = sorted(list(data['events']['tab-size'].unique()))

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
    for ts in ts_keys:
        for t, topology in enumerate(topology_keys):

            _data = data['events'].loc[
                (data['events']['topology'] == topology)
                & (data['events']['tab-size'] == ts)]

            label = ['', '']
            if setlabel == True:
                label[0] = topology_labels[topology]
                if t == 0:
                    label[1] = 'fwd. effic.'

            ax2.bar(
                pos, 
                _data['llm'] + _data['slm'] - _data['mlm'],
                color = topology_colors[t], linewidth = 0.5, width = bar_width, label = label[0])

            ax1.plot(
                pos + (bar_width / 2.0),
                ((4.0 / (_data['llm'] + _data['slm'] - _data['mlm'])) * 2.0),
                linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

            if t == 2:
                xbf[pos] = '10^' + str(int(math.log10(ts)))

            pos += bar_width

        # only include lable once
        setlabel = False
        pos += (bar_width * 2.0)

    # an horizontal shaded area, for fwd efficiency
    ax1.axhspan(0, 2, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    # legend
    ax2.legend(fontsize = 10, ncol = 1, loc = 'upper right')
    ax1.legend(fontsize = 10, ncol = 1, loc = 'upper left')

    # axis labels
    ax1.set_xlabel("fwd. table size\n\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax2.set_ylabel("avg. # of used links")
    ax1.set_ylabel("fwd. efficiency [0,1]")

    # avg used link axis
    ax2.set_xlim(0.0 - bar_width, pos - (1.0 * bar_width))
    ax2.set_ylim(1, 10000000)
    ax2.set_yscale("log", nonposy = 'clip')
    ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    
    # fwd effic axix
    ax1.set_ylim(-2, 5)
    ax1.set_yticks([0, 1, 2])
    ax1.set_yticklabels([0, 0.5, 1])

    plt.savefig(os.path.join(output_dir, "global-tab-sizes.pdf"), bbox_inches='tight', format = 'pdf')

def plot_bf_sizes(data_dir, test_dir, output_dir):

    # get data
    data = read_data(data_dir, ['events'])
    for file_type in ['events']:
        print(data[file_type])

    # labels & stuff
    topology_keys = sorted(list(data['events']['topology'].unique()))
    topology_keys.append('ideal')
    # grayscale for topologies
    topology_colors = ['#000000', '#bebebe', 'white', '#708090']
    bf_keys = sorted(list(data['events']['bf-size'].unique()))
    fwd_entry_keys = sorted(list(data['events']['fwd-size'].unique()))

    print(fwd_entry_keys)

    # prepare the graph parameters
    bar_width = 0.10

    # use the classic plot style
    plt.style.use('classic')
    # avoid Type 3 fonts
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # fig
    fig = plt.figure(figsize=(5, 3.5))
    # (1) avg. nr. of used links (bar chart)
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True, which = 'major')

    # (2) fwd efficiency (line chart)
    # ax2 = ax1.twinx()
    # ax2.yaxis.grid(True, ls = 'dotted', lw = 0.50)

    # ax1.set_zorder(ax2.get_zorder() + 1)
    # ax1.patch.set_visible(False)

    pos = 0.0
    setlabel = True
    xbf = defaultdict()
    for t, topology in enumerate(topology_keys):

        if topology == 'ideal':
            ax1.axvline(x = (pos - (0.75 * bar_width)), linewidth = 0.5, linestyle = '-', color = 'black')
            ax1.axvspan(pos - (0.75 * bar_width), pos + (8 * bar_width) + bar_width, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.50)

        for f, fwd in enumerate(fwd_entry_keys):
            for b, bf in enumerate(bf_keys):

                value = 0.0
                if topology != 'ideal':
                    _data = data['events'].loc[
                        (data['events']['topology'] == topology)
                        & (data['events']['bf-size'] == bf)
                        & (data['events']['fwd-size'] == fwd)]

                    value = _data['mlm'] + _data['slmc'] + _data['slmi']

                else:
                    value = 4.0

                label = ['', '']
                if setlabel == True:
                    label[0] = bf
                    if b == 0:
                        label[1] = 'fwd. effic.'

                ax1.bar(
                    pos, 
                    value,
                    color = topology_colors[b], linewidth = 0.5, width = bar_width, label = label[0])

                # ax1.plot(
                #     pos + (bar_width / 2.0),
                #     ((4.0 / (_data['llm'] + _data['slm'] - _data['mlm'])) * 2.0),
                #     linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

                if b == 1:
                    xbf[pos + (bar_width)] = str(fwd)

                pos += bar_width

            if f == 0:
                xbf[pos + (0.5 * bar_width)] = '\n' + str(topology_labels[topology])

            # only include lable once
            setlabel = False
            pos += (bar_width / 2.0)

        pos += (bar_width * 1.0)

    # an horizontal shaded area, for fwd efficiency

    # legend
    leg = []
    leg.append(ax1.legend(
        fontsize = 12, 
        ncol = 2, loc = 'upper right', title = 'bf size',
        handletextpad = 0.2, handlelength = 1.0, labelspacing = 0.2, columnspacing = 0.5))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 12)

    # ax1.legend(fontsize = 10, ncol = 1, loc = 'upper left')

    # axis labels
    ax1.set_xlabel("fwd. entry size\ntopology")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax1.set_ylabel("avg. # of used links")
    # ax1.set_ylabel("fwd. efficiency [0,1]")

    # avg used link axis
    ax1.set_xlim(0.0 - (bar_width / 2.0), pos - (1.0 * bar_width))
    ax1.set_ylim(1, 10000)
    ax1.set_yscale("log", nonposy = 'clip')
    ax1.set_yticks([10**(x) for x in np.arange(0, 5, 1)])
    ax1.tick_params(axis = 'y', which = 'minor', bottom = 'off')
    
    # # fwd effic axix
    # ax1.set_ylim(-2, 5)
    # ax1.set_yticks([0, 1, 2])
    # ax1.set_yticklabels([0, 0.5, 1])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "global-bf-sizes.pdf"), format = 'pdf') 

