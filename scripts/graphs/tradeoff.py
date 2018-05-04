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
            __data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst', 'type', 'latency'])['prob'].agg('sum').reset_index()
            __data['avg-latency'] = __data['latency'] * __data['prob']

            # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size', 'src:dst', 'type'}, and sum() by 'prob'
            _data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst', 'type'])['prob'].agg('sum').to_frame().reset_index()
            _data['avg-latency'] = __data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst', 'type'])['avg-latency'].agg('sum').reset_index()['avg-latency'] / _data['prob']

        elif file_type == 'events':
            # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size', 'src:dst'}, and sum() the diff event types
            _data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst'])['llm', 'mlm', 'slm'].agg('sum').reset_index()

        # update main dataframe
        data[file_type] = pd.concat([data[file_type], _data], ignore_index = True)

    # calc sum of each diff. outcome (we use the term 'event' in the paper),
    # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size'} 
    for file_type in file_types:
        if file_type == 'outcomes':
            data[file_type] = data[file_type].groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'type'])[['prob', 'avg-latency']].agg('mean').reset_index()
        elif file_type == 'events':
            data[file_type] = data[file_type].groupby(['topology', 'bf-size', 'req-size', 'tab-size', 'fwd-size'])[['llm', 'mlm', 'slm']].agg('mean').reset_index()

    return data

def plot_tradeoff(data_dir, test_dir, output_dir):
    
    data = read_data(data_dir, ['events', 'outcomes'])

    bf_keys = sorted(list(data['events']['bf-size'].unique()))
    rs_keys = sorted(list(data['events']['req-size'].unique()))

    print(bf_keys)

    linestyles = ['-','--']
    markers = ['o', '^', 'v', 's', 'x']

    # use the classic plot style
    plt.style.use('classic')

    # fig
    fig = plt.figure(figsize=(5, 3.5))
    xlims = [0.5, 0.5]
    ylims = [[0.0, 1.0], [0.0, 1.0]]
    ylabels = ['bf fwd. use (%)', 'fallback use (%)']

    ax = defaultdict()
    _data = defaultdict()
    for rr, rs in enumerate(rs_keys):
        
        _data['events'] = data['events'].loc[(data['events']['req-size'] == rs)]
        _data['outcomes'] = data['outcomes'].loc[(data['outcomes']['req-size'] == rs)]
        
        series = defaultdict(list)

        for bb, bf in enumerate(bf_keys):

            _d  = _data['events'].loc[_data['events']['bf-size'] == bf]
            __d = _data['outcomes'].loc[(_data['outcomes']['bf-size'] == bf) & (_data['outcomes']['type'] == 'hfd')]

            print("%s, %s" % (rs, bf))
            print("MLM : %s" % (_d['mlm'].values))
            print("SLM : %s" % (_d['slm'].values))

            hfd_value = 0.0
            if not __d.empty:
                print(__d)
                hfd_value = __d['prob'].values * (4.0 - __d['avg-latency'].values)
                print("HFD : %s" % (hfd_value))
                print("HFD (lat) : %s" % (__d['avg-latency'].values))

            total = _d['mlm'].values + _d['slm'].values + hfd_value
            print("TOTAL : %s" % (total))
            series['bf'].append((_d['mlm'].values + _d['slm'].values) / total)
            series['fb'].append(hfd_value / total)

        for c, case in enumerate(['bf', 'fb']):

            ax[c] = fig.add_subplot(211 + c)
            ax[c].xaxis.grid(True)
            ax[c].yaxis.grid(True)

            ax[c].plot(
                np.arange(0, len(series[case]), 1),
                series[case],
                linewidth = 1.0, color = 'black', linestyle = '-', marker = markers[rr], markersize = 5, 
                label = rs)

            ax[c].set_ylabel(ylabels[c])
            ax[c].set_ylim(ylims[c][0], ylims[c][1])
            ax[c].set_yticks(np.arange(ylims[c][0], ylims[c][1] + 0.2, 0.2))

            ax[c].set_xlim(0 - 0.1, (len(bf_keys) - 1) + 0.1)
            ax[c].set_xticks(np.arange(0, len(bf_keys), 1))
            ax[c].set_xticklabels(bf_keys)

    # legend
    leg = []
    leg.append(ax[0].legend(fontsize = 10, ncol = 2, loc = 'lower right', title = 'req. size', handletextpad = 0.2))
    leg.append(ax[1].legend(fontsize = 10, ncol = 2, loc = 'upper right', title = 'req. size', handletextpad = 0.2))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 10)

    # title
    ax[0].set_title("AT&T (7018)", fontsize = 10)
    # axis labels
    ax[1].set_xlabel("bf size")

    # fig.subplots_adjust(left = None, bottom = None, right = None, top = None, wspace = None, hspace = 0.40)
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, "global-tradeoff.pdf"), format = 'pdf')
