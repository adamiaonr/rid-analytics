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
            _data = _data.groupby(['topology', 'bf-size', 'req-size', 'fwd-size', 'tab-size', 'src:dst'])['llm', 'mlm', 'slmc', 'slmi', 'fpd'].agg('sum').reset_index()

        # update main dataframe
        data[file_type] = pd.concat([data[file_type], _data], ignore_index = True)

    # calc sum of each diff. outcome (we use the term 'event' in the paper),
    # groupby() {'topology', 'bf-size', 'req-size', 'fwd-size'} 
    for file_type in file_types:
        if file_type == 'outcomes':
            data[file_type] = data[file_type].groupby(['type'])[['prob', 'avg-latency']].agg('mean').reset_index()
        elif file_type == 'events':
            data[file_type] = data[file_type].groupby(['topology', 'bf-size', 'req-size', 'tab-size', 'fwd-size'])[['llm', 'mlm', 'slmc', 'slmi', 'fpd']].agg('mean').reset_index()

    return data

def plot_breakdown(data_dir, test_dir, scenario, output_dir):
    
    # extract diff data dirs
    data_dir = data_dir.split(':')
    # grayscale colors
    colors = ['lightgray', 'black', '#708090']

    # labels, keys & other stuff
    topology_keys = ['7018']
    # grayscale colors
    colors = ['lightgray', '#000000', '#708090', 'white']
    # mm resolution modes    
    modes = [0, 3]
    modes_str = {'0192' : 'AML\n(192)', '0384' : 'AML (384 bit)', '3192' : 'FB (192 bit)'}
    # bf sizes
    bf_sizes = ['192', '384']

    cases = []
    cstr = defaultdict()
    if scenario == 'cdn':
        cases = ['S510021', 'S510023', 'S51002I']
    elif scenario == 'opportunistic':
        cases = ['S1010021', 'S1010023', 'S101002I']
    elif scenario == 'cdn-lpm':
        cases = ['S350221S450221', 'S35022IS45022I', 'S450221S550221', 'S45022IS55022I', 'S550221', 'S55022I']

    # prepare the graph parameters
    bar_width = 0.10
    # use the classic plot style
    plt.style.use('classic')
    # avoid Type 3 fonts
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # fig
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(121)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True, ls = 'dotted', lw = 0.50)

    # aux. variables
    pos = 0.0
    xbf = defaultdict()
    ### 1st group : caching scenario (no tps)
    ax1.set_title("(a) no TPs", fontsize = 10)
    for t, topology in enumerate(topology_keys):
        for case in cases:

            xx = (bar_width / 2.0)
            for bf in bf_sizes:

                print("topology : %s, case : %s, bf-size : %s" % (topology, case, bf))
                _data = read_data(data_dir[0], file_types = ['events'], filters = {'topology' : topology, 'case' : case, 'bf-size' : bf})['events']

                # fwd corr. cases
                fwd_corr = defaultdict(float)
                # fwd corr. labels
                labels = {'slmc' : '', 'slmi' : '', 'fpd' : ''}

                # sum of fwd corr. case values for % calculation
                total = _data['slmi'].values + _data['slmc'].values + _data['fpd'].values
                prev = 0.0
                for fc, fwd_case in enumerate(['slmc', 'slmi', 'fpd']):
                    fwd_corr[fwd_case] = (_data[fwd_case].values / total) * 100.0
                    print('%s : %s' % (fwd_case, fwd_corr[fwd_case]))

                    # print stacked bar
                    _bar = ax1.bar(
                        pos, 
                        fwd_corr[fwd_case],
                        color = colors[fc], linewidth = 0.50, width = bar_width, label = labels[fwd_case],
                        bottom = prev)

                    # offset for next stacked portion of the bar
                    prev += fwd_corr[fwd_case]

                    # print bf size

                    if (t == 0) and (case == cases[0]) and (fwd_case == 'fpd'):
                        for rect in _bar:
                            height = prev
                            plt.text(rect.get_x() + xx, height, str(bf) + ' bit', ha = 'center', va = 'bottom', rotation = 80, fontsize = 10)
                            xx += (rect.get_width())

                # print the # neighboring caches case
                if (bf == '192'):
                    # xbf[pos + (bar_width)] = '\n' + str(case[1] + ':' + (case[-1] if (case[-1] != 'I') else 'all'))
                    xbf[pos + (bar_width)] = str((case[-1] if (case[-1] != 'I') else '*'))

                pos += bar_width
            pos += (bar_width / 2.0)

            if (case == cases[1]):
                # xbf[pos + (bar_width)] = '\n' + str(case[1] + ':' + (case[-1] if (case[-1] != 'I') else 'all'))
                xbf[pos - (1.499 * bar_width)] = '\n' + str(topology_labels[int(topology)])

        pos += (bar_width * 1.0)

    # axis labels
    ax1.set_xlabel("\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())
    # for tick in ax1.get_xticklabels():
    #     tick.set_rotation(45)
    ax1.set_ylabel("share of fwd.\n decision type (%)")
    # ax1.set_ylabel("fwd. efficiency [0,1]")
    # avg used link axis
    ax1.set_xlim(0.0 - (bar_width), pos - (0.5 * bar_width))
    # fwd effic axix
    ax1.set_ylim(0.0, 140.0)
    ax1.set_yticks(np.arange(0.0, 110.0, 20.0))

    ### 2nd group : caching scenario (w/ tps)

    ax1 = fig.add_subplot(122)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True, ls = 'dotted', lw = 0.50)
    ax1.set_title("(b) w/ TPs", fontsize = 10)
    # aux. variables
    pos = 0.0
    xbf = defaultdict()
    ### 2nd group : caching scenario (w/ tps)

    if scenario == 'cdn-lpm':
        cases = ['S35022IS45022I*132', 'S35022IS45022I*142', 'S45022IS55022I*142', 'S45022IS55022I*152']

    for t, topology in enumerate(topology_keys):
        for case in cases:
            for bf in bf_sizes:

                print("topology : %s, case : %s, bf-size : %s" % (topology, case, bf))
                _data = read_data(data_dir[1], file_types = ['events'], filters = {'topology' : topology, 'case' : case, 'bf-size' : bf})['events']

                # fwd corr. cases
                fwd_corr = defaultdict(float)
                # fwd corr. labels 
                if (t == 0) and (case == cases[0]) and (bf == '192'):
                    labels = {'slmc' : 'corr.', 'slmi' : 'incorr.', 'fpd' : 'fp detected'}
                else:
                    labels = {'slmc' : '', 'slmi' : '', 'fpd' : ''}

                # sum of fwd corr. case values for % calculation
                total = _data['slmi'].values + _data['slmc'].values + _data['fpd'].values
                prev = 0.0
                for fc, fwd_case in enumerate(['slmc', 'slmi', 'fpd']):
                    fwd_corr[fwd_case] = (_data[fwd_case].values / total) * 100.0
                    print('%s : %s' % (fwd_case, fwd_corr[fwd_case]))

                    # print stacked bar
                    _bar = ax1.bar(
                        pos, 
                        fwd_corr[fwd_case],
                        color = colors[fc], linewidth = 0.50, width = bar_width, label = labels[fwd_case],
                        bottom = prev)

                    # offset for next stacked portion of the bar
                    prev += fwd_corr[fwd_case]

                # print the # neighboring caches case
                if (bf == '192'):
                    # xbf[pos + (bar_width)] = '\n' + str(case[1] + ':' + (case[-1] if (case[-1] != 'I') else 'all'))
                    xbf[pos + (bar_width)] = str((case[-1] if (case[-1] != 'I') else '*'))

                pos += bar_width
            pos += (bar_width / 2.0)

            if (case == cases[1]):
                # xbf[pos + (bar_width)] = '\n' + str(case[1] + ':' + (case[-1] if (case[-1] != 'I') else 'all'))
                xbf[pos - (1.499 * bar_width)] = '\n' + str(topology_labels[int(topology)])

        pos += (bar_width * 1.0)

    # legend
    leg = []
    leg.append(ax1.legend(fontsize = 10, ncol = 3, loc = 'upper right', title = 'fwd. decision correctness', handletextpad = 0.2, bbox_to_anchor=(0.900, 1.0)))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 10)

    # axis labels
    # ax1.set_xlabel("# of neighboring caches\ntopology")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    fig.text(0.575, 0.05, "# of caches around req. source\ntopology", ha = 'center')

    # avg used link axis
    ax1.set_xlim(0.0 - (bar_width), pos - (0.5 * bar_width))
    
    # fwd effic axix
    ax1.set_ylim(0.0, 140.0)
    ax1.set_yticks(np.arange(0.0, 110.0, 20.0))
    ax1.set_yticklabels([''] * 6)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, ("%s-breakdown.pdf" % (scenario))), format = 'pdf')

def plot_usage(data_dir, test_dir, scenario, output_dir):

    # extract diff data dirs
    data_dir = data_dir.split(':')

    # labels, keys & other stuff
    topology_keys = ['7018', '3257']
    # grayscale colors
    colors = ['#000000', '#708090', '#bebebe', 'lightgray']
    # mm resolution modes    
    modes = [0, 3]
    modes_str = {'0192' : 'AML\n(192)', '0384' : 'AML (384 bit)', '3192' : 'FB (192 bit)'}
    # bf sizes (per mm resolution mode)
    bf_sizes = {0 : 384, 3 : 192}

    cases = []
    cstr = defaultdict()
    if scenario == 'cdn':
        cases = ['S212', 'S512']
        cstr = {'S212' : 2, 'S512' : 5}

    elif scenario == 'opportunistic':
        cases = ['S101002I']
        cstr = {'S101002I' : 10}

    # graph parameters
    bar_width = 0.10
    # use the classic plot style
    plt.style.use('classic')
    # avoid Type 3 fonts
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # fig
    fig = plt.figure(figsize=(5, 3.5))
    # (1) avg. nr. of used links (bar chart)
    ax1 = fig.add_subplot(121)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(False)
    # (2) latency (line chart)
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, ls = 'dotted', lw = 0.50)
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    pos = 0.0
    setlabel = True
    xbf = defaultdict()
    ax1.set_title("(a) no TPs", fontsize = 10)
    for t, topology in enumerate(topology_keys):
        for c in cases:
            for m in modes:

                bf = bf_sizes[m]
                print("topology : %s, mode : %s, case : %s, bf-size : %s" % (topology, m, c, bf))

                cd_lat = 0.0
                link_usage = 0.0
                if topology != '1337':

                    _data = read_data(data_dir[0], file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : c, 'mode' : m})

                    # bars  : link usage
                    # lines : avg. latency
                    cd_lat = float(_data['outcomes'].loc[_data['outcomes']['type'] == 'cd']['avg-latency'])

                    # SLM events are invariable for either mode
                    link_usage = _data['events']['slmi'].values + _data['events']['slmc'].values

                    # regarding MLM events:
                    #   - if mode is AML (m = 0), then we simply sum the avg. nr. of mlm events
                    #   - if mode is FB (m = 3), we must some the LLM events + the fwd events due 
                    #     to fallbacks, which are approx. the avg. latency of fallback events
                    if m == 0:
                        link_usage = (link_usage + _data['events']['mlm'].values) / 4.0
                    if m == 3:
                        __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']

                        print("link_usage = %s + %s + %s = %s" % 
                            (link_usage, _data['events']['llm'].values, 
                            __d['avg-latency'].values, 
                            link_usage + _data['events']['llm'].values + (__d['avg-latency'].values)))
                        link_usage = (link_usage + (__d['prob'].values * __d['avg-latency'].values)) / 4.0

                else:
                    cd_lat = 4.0
                    link_usage = (4.0 / 4.0)

                label = ['', '', '']
                if setlabel == True:
                    label[0] = str(modes_str[str(m) + str(bf)])
                    if m == 0:
                        label[1] = 'avg. latency'

                ax2.bar(
                    pos, 
                    link_usage,
                    color = colors[m], linewidth = 0.5, width = bar_width, label = label[0])

                # lines hold avg. # of deliveries
                ax1.plot(
                    pos + (bar_width / 2.0),
                    cd_lat,
                    linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

                if m == 0:
                    xbf[pos + (bar_width)] = str(cstr[c])

                pos += bar_width

            if c == cases[0]:
                if c == 'S212':
                    xbf[pos] = '\n' + str(topology_labels[int(topology)])
                elif c == 'S1010021': 
                    xbf[pos - (0.99 * bar_width)] = '\n' + str(topology_labels[int(topology)])

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
        [4.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    ax1.plot(
        [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        [5.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    # legends
    leg = []
    leg.append(ax1.legend(fontsize = 10, ncol = 1, loc = 'upper left'))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 10)

    # axis labels
    ax1.set_xlabel("\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    # ax2.set_ylabel("avg. # of used links")
    ax1.set_ylabel("avg. latency")

    # avg used link axis
    ax2.set_xlim(0.0 - (bar_width / 2.0), pos - (1.0 * bar_width))
    ax2.set_ylim(0, 5)
    ax2.set_yticks([0, 1, 2])
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    # ax2.set_ylim(1, 10000)
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 5, 1)])
    # ax2.tick_params(axis = 'y', which = 'minor', bottom = 'off')
    
    # fwd effic axis
    ax1.set_ylim(-2, 8)
    ax1.set_yticks([3.0, 4.0, 5.0])

    # (1) avg. nr. of used links (bar chart)
    ax1 = fig.add_subplot(122)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(False)
    # (2) latency (line chart)
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, ls = 'dotted', lw = 0.50)
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    pos = 0.0
    setlabel = True
    xbf = defaultdict()
    ax1.set_title("(b) w/ TPs", fontsize = 10)
    for t, topology in enumerate(topology_keys):
        for c in cases:
            for m in modes:

                bf = bf_sizes[m]
                print("topology : %s, mode : %s, case : %s, bf-size : %s" % (topology, m, c, bf))

                cd_lat = 0.0
                link_usage = 0.0
                if topology != '1337':

                    _data = read_data(data_dir[1], file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : c, 'mode' : m})

                    # bars  : link usage
                    # lines : avg. latency
                    cd_lat = float(_data['outcomes'].loc[_data['outcomes']['type'] == 'cd']['avg-latency'])

                    # SLM events are invariable for either mode
                    link_usage = _data['events']['slmi'].values + _data['events']['slmc'].values

                    # regarding MLM events:
                    #   - if mode is AML (m = 0), then we simply sum the avg. nr. of mlm events
                    #   - if mode is FB (m = 3), we must some the LLM events + the fwd events due 
                    #     to fallbacks, which are approx. the avg. latency of fallback events
                    if m == 0:
                        link_usage = (link_usage + _data['events']['mlm'].values) / 2.0
                    if m == 3:
                        __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']

                        print("link_usage = %s + %s + %s = %s" % 
                            (link_usage, _data['events']['llm'].values, 
                            __d['avg-latency'].values, 
                            link_usage + _data['events']['llm'].values + (__d['avg-latency'].values)))
                        link_usage = (link_usage + (__d['prob'].values * __d['avg-latency'].values)) / 2.0

                else:
                    cd_lat = 2.0
                    link_usage = (2.0 / 2.0)

                label = ['', '', '']
                if setlabel == True:
                    label[0] = str(modes_str[str(m) + str(bf)])
                    if m == 0:
                        label[1] = 'avg. latency'

                ax2.bar(
                    pos, 
                    link_usage,
                    color = colors[m], linewidth = 0.5, width = bar_width, label = label[0])

                # lines hold avg. # of deliveries
                ax1.plot(
                    pos + (bar_width / 2.0),
                    cd_lat,
                    linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = label[1])

                if m == 0:
                    xbf[pos + (bar_width)] = str(cstr[c])

                pos += bar_width

            if c == cases[0]:
                if c == 'S212':
                    xbf[pos] = '\n' + str(topology_labels[int(topology)])
                elif c == 'S1010021': 
                    xbf[pos - (0.99 * bar_width)] = '\n' + str(topology_labels[int(topology)])

            # only include lable once
            setlabel = False
            pos += (bar_width / 2.0)

        pos += (bar_width * 1.0)

    # an horizontal shaded area, for fwd efficiency
    ax1.axhspan(2, 4, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.25)

    # ax1.plot(
    #     [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
    #     [1.0] * 2,
    #     linewidth = 0.5, color = 'black', linestyle = ':')

    ax1.plot(
        [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        [2.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    ax1.plot(
        [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        [3.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    ax1.plot(
        [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        [4.0] * 2,
        linewidth = 0.5, color = 'black', linestyle = ':')

    # legends
    leg = []
    leg.append(ax2.legend(fontsize = 10, ncol = 1, loc = 'upper right', title = 'fwd. strategy'))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 10)

    # axis labels
    ax1.set_xlabel("\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax2.set_ylabel("link usage excess")
    # ax1.set_ylabel("avg. latency")

    fig.text(0.50, 0.05, "cache annc. size\ntopology", ha = 'center')

    # avg used link axis
    ax2.set_xlim(0.0 - (bar_width / 2.0), pos - (1.0 * bar_width))
    ax2.set_ylim(0, 5)
    ax2.set_yticks([0, 1, 2, 3])
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    # ax2.set_ylim(1, 10000)
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 5, 1)])
    # ax2.tick_params(axis = 'y', which = 'minor', bottom = 'off')
    
    # fwd effic axis
    ax1.set_ylim(-2, 8)
    ax1.set_yticks([2.0, 3.0, 4.0])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, ("%s-usage.pdf" % (scenario))), bbox_inches='tight', format = 'pdf')

def plot_breakdown_(data_dir, test_dir, scenario, output_dir):
    
    # extract diff data dirs
    data_dir = data_dir.split(':')
    # grayscale colors
    colors = ['lightgray', 'black', '#708090']

    # labels, keys & other stuff
    topology_keys = ['3356', '7018']
    # grayscale colors
    colors = ['lightgray', 'dimgray', 'black', 'white']
    # mm resolution modes    
    modes = [0, 3]
    modes_str = {'384' : 'AML', '192' : 'FB'}
    # bf sizes
    bf_sizes = ['384', '192']

    cases = OrderedDict()
    cstr = defaultdict()

    if scenario == 'opportunistic':

        metacases = ['{7,10}']
        _cases = {
            # '{7,10}' : ['S750221S1050221', 'S75022HS105022H', 'S75022IS105022I']
            '{7,10}' : ['S75022HS105022H', 'S75022IS105022I']
            }

    elif scenario == 'cdn-lpm':

        metacases = ['{4,5}']
        _cases = {
            '{2,3}' : ['S250221S350221', 'S25022HS35022H', 'S25022IS35022I'], 
            '{3,4}' : ['S350221S450221', 'S35022HS45022H', 'S35022IS45022I'], 
            # '{4,5}' : ['S450221S550221', 'S45022HS55022H', 'S45022IS55022I'] 
            '{4,5}' : ['S45022HS55022H', 'S45022IS55022I'] 
            }

    for metacase in metacases:
        cases[metacase] = _cases[metacase]

    # prepare the graph parameters
    bar_width = 0.10
    # use the classic plot style
    plt.style.use('classic')
    # avoid Type 3 fonts
    matplotlib.rcParams['pdf.fonttype'] = 42
    matplotlib.rcParams['ps.fonttype'] = 42

    # fig
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(121)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True, ls = 'dotted', lw = 0.50)

    # aux. variables
    pos = 0.0
    xbf = defaultdict()
    ### 1st group : caching scenario (no tps)
    ax1.set_title("(a) no TPs", fontsize = 12)
    for t, topology in enumerate(topology_keys):
        for bf in bf_sizes:
            for metacase in cases:
                for case in cases[metacase]:

                    print("topology : %s, case : %s, bf-size : %s" % (topology, case, bf))
                    _data = read_data(data_dir[0], file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : case, 'bf-size' : bf})

                    # fwd corr. cases
                    fwd_corr = defaultdict(float)
                    # fwd corr. labels
                    labels = {'slmc' : '', 'slmi' : '', 'fpd' : ''}

                    hfd_contrib = 0.0
                    # if bf == '192':
                    #     __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']
                    #     hfd_contrib = (__d['prob'].values * __d['avg-latency'].values)

                    # sum of fwd corr. case values for % calculation
                    total = _data['events']['slmi'].values + _data['events']['slmc'].values + _data['events']['fpd'].values + hfd_contrib
                    prev = 0.0
                    for fc, fwd_case in enumerate(['slmc', 'slmi', 'fpd']):
                        if fwd_case == 'slmc':
                            fwd_corr[fwd_case] = ((_data['events'][fwd_case].values + hfd_contrib) / total) * 100.0
                        else:
                            fwd_corr[fwd_case] = (_data['events'][fwd_case].values / total) * 100.0
                        print('%s : %s' % (fwd_case, fwd_corr[fwd_case]))

                        # print stacked bar
                        _bar = ax1.bar(
                            pos, 
                            fwd_corr[fwd_case],
                            color = colors[fc], linewidth = 0.50, width = bar_width, label = labels[fwd_case],
                            bottom = prev)

                        # offset for next stacked portion of the bar
                        prev += fwd_corr[fwd_case]

                        # print # of caches case
                        if (t == 0) and (bf == '384') and (metacase == cases.keys()[0]) and (fwd_case == 'fpd'):
                            for rect in _bar:
                                height = prev + 1
                                plt.text(
                                    rect.get_x() + (bar_width / 2.0), height, 
                                    ('n' if (case[-1] == 'I') else (('n/2 caches') if (case[-1] == 'H') else ('%s cache' % (case[-1])))), 
                                    ha = 'center', va = 'bottom', rotation = 75, fontsize = 10)
                    
                    pos += bar_width
                    if (case == cases[metacase][0]):
                        xbf[pos] = str(modes_str[bf])

                if (bf == '384'):
                    xbf[pos + (bar_width / 2.0)] = '\n' + str(topology_labels[int(topology)])

                pos += (bar_width / 2.0)
            pos += (bar_width / 2.0)
        pos += (bar_width * 1.0)

    # axis labels
    ax1.set_xlabel("\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())
    # for tick in ax1.get_xticklabels():
    #     tick.set_rotation(45)
    ax1.set_ylabel("% of fwd. decisions,\nby correctness")
    # ax1.set_ylabel("fwd. efficiency [0,1]")
    # avg used link axis
    ax1.set_xlim(0.0 - (bar_width), pos - (1.0 * bar_width))
    # fwd effic axix
    ax1.set_ylim(0.0, 150)
    ax1.set_yticks(np.arange(0.0, 110.0, 20.0))

    ### 2nd group : caching scenario (w/ tps)

    ax1 = fig.add_subplot(122)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True, ls = 'dotted', lw = 0.50)
    ax1.set_title("(b) w/ TPs", fontsize = 12)
    # aux. variables
    pos = 0.0
    xbf = defaultdict()
    ### 2nd group : caching scenario (w/ tps)

    cases = OrderedDict()
    cstr = defaultdict()
    if scenario == 'opportunistic':

        metacases = ['{7,10}']
        _cases = {
            # '{7,10}' : ['S750221S1050221*172', 'S75022HS105022H*172', 'S75022IS105022I*172']
            '{7,10}' : ['S75022HS105022H*172', 'S75022IS105022I*172']
            }

    elif scenario == 'cdn-lpm':

        metacases = ['{4,5}']
        _cases = {
            '{2,3}' : ['S250221S350221*132', 'S25022HS35022H*132', 'S25022IS35022I*132'], 
            '{3,4}' : ['S350221S450221*142', 'S35022HS45022H*142', 'S35022IS45022I*142'], 
            # '{4,5}' : ['S450221S550221*152', 'S45022HS55022H*152', 'S45022IS55022I*152']
            '{4,5}' : ['S45022HS55022H*152', 'S45022IS55022I*152']
            }

    for metacase in metacases:
        cases[metacase] = _cases[metacase]

    for t, topology in enumerate(topology_keys):
        for bf in bf_sizes:
            for metacase in cases:
                for case in cases[metacase]:

                    print("topology : %s, case : %s, bf-size : %s" % (topology, case, bf))
                    _data = read_data(data_dir[1], file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : case, 'bf-size' : bf})

                    # fwd corr. cases
                    fwd_corr = defaultdict(float)
                    # fwd corr. labels
                    labels = None
                    if (t == 0) and (metacase == cases.keys()[0]) and (case == cases[metacase][0]) and (bf == '192'):
                        labels = {'slmc' : 'corr.', 'slmi' : 'incorr.', 'fpd' : 'mult. match'}
                    else:
                        labels = {'slmc' : '', 'slmi' : '', 'fpd' : ''}

                    hfd_contrib = 0.0
                    # if bf == '192':
                    #     __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']
                    #     hfd_contrib = (__d['prob'].values * __d['avg-latency'].values)

                    # sum of fwd corr. case values for % calculation
                    total = _data['events']['slmi'].values + _data['events']['slmc'].values + _data['events']['fpd'].values + hfd_contrib
                    prev = 0.0
                    for fc, fwd_case in enumerate(['slmc', 'slmi', 'fpd']):
                        if fwd_case == 'slmc':
                            fwd_corr[fwd_case] = ((_data['events'][fwd_case].values + hfd_contrib) / total) * 100.0
                        else:
                            fwd_corr[fwd_case] = (_data['events'][fwd_case].values / total) * 100.0
                        print('%s : %s' % (fwd_case, fwd_corr[fwd_case]))

                        # print stacked bar
                        _bar = ax1.bar(
                            pos, 
                            fwd_corr[fwd_case],
                            color = colors[fc], linewidth = 0.50, width = bar_width, label = labels[fwd_case],
                            bottom = prev)

                        # offset for next stacked portion of the bar
                        prev += fwd_corr[fwd_case]

                        # print # of caches case
                        # if (t == 0) and (bf == '384') and (metacase == cases.keys()[0]) and (fwd_case == 'fpd'):
                        #     for rect in _bar:
                        #         height = prev + 1
                        #         plt.text(
                        #             rect.get_x() + (bar_width / 2.0), height, 
                        #             ('n' if (case[-1] == 'I') else (('n/2') if (case[-1] == 'H') else ('%s cache' % (case[-1])))), 
                        #             ha = 'center', va = 'bottom', rotation = 75, fontsize = 10)
                    
                    pos += bar_width
                    if (case == cases[metacase][0]):
                        xbf[pos] = str(modes_str[bf])

                if (bf == '384'):
                    xbf[pos + (bar_width / 2.0)] = '\n' + str(topology_labels[int(topology)])

                pos += (bar_width / 2.0)
            pos += (bar_width / 2.0)
        pos += (bar_width * 1.0)

    # legend
    leg = []
    leg.append(ax1.legend(
        fontsize = 12, 
        ncol = 3, loc = 'upper right', 
        title = 'fwd. decision correctness', 
        handletextpad = 0.2, handlelength = 1.25, labelspacing = 0.2, columnspacing = 0.5,
        bbox_to_anchor=(0.825, 1.000)))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 12)

    # axis labels
    # ax1.set_xlabel("# of neighboring caches\ntopology")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    fig.text(0.575, 0.065, "fwd. strategy (bf size)\ntopology", ha = 'center')

    # avg used link axis
    ax1.set_xlim(0.0 - (bar_width), pos - (1.0 * bar_width))
    
    # fwd effic axix
    ax1.set_ylim(0.0, 150)
    ax1.set_yticks(np.arange(0.0, 110.0, 20.0))
    ax1.set_yticklabels([''] * 6)

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, ("%s-breakdown.pdf" % (scenario))), format = 'pdf') 

def plot_usage_(data_dir, test_dir, scenario, output_dir):

    # extract diff data dirs
    data_dir = data_dir.split(':')

    # labels, keys & other stuff
    topology_keys = ['3356', '7018']
    # grayscale colors
    colors = ['#000000', '#708090', '#bebebe', 'lightgray']
    # mm resolution modes    
    modes = [0, 3]
    modes_str = {'384' : 'AML(384)', '192' : 'FB(192)'}
    # bf sizes
    bf_sizes = {3 : '192', 0 : '384'}

    markers = {'384' : '^', '192' : 'o'}

    cases = OrderedDict()
    cstr = defaultdict()

    if scenario == 'opportunistic':

        metacases = ['{7,10}']
        _cases = {
            # '{7,10}' : ['S750221S1050221', 'S75022HS105022H', 'S75022IS105022I']
            '{7,10}' : ['S75022HS105022H', 'S75022IS105022I']
            }

    elif scenario == 'cdn-lpm':

        metacases = ['{4,5}']
        _cases = {
            '{2,3}' : ['S250221S350221', 'S25022HS35022H', 'S25022IS35022I'], 
            '{3,4}' : ['S350221S450221', 'S35022HS45022H', 'S35022IS45022I'], 
            # '{4,5}' : ['S450221S550221', 'S45022HS55022H', 'S45022IS55022I'] 
            '{4,5}' : ['S45022HS55022H', 'S45022IS55022I'] 
            }

    for metacase in metacases:
        cases[metacase] = _cases[metacase]

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
    ax1 = fig.add_subplot(121)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(False)
    # (2) latency (line chart)
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, ls = 'dotted', lw = 0.50)
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    pos = 0.0
    setlabel = True
    xbf = defaultdict()
    ax1.set_title("(a) no TPs", fontsize = 12)

    latency = defaultdict()

    for t, topology in enumerate(topology_keys):

        latency[0] = defaultdict(list)
        latency[3] = defaultdict(list)

        for metacase in cases:
            for case in cases[metacase]:
                for m in modes:

                    bf = bf_sizes[m]
                    print("topology : %s, mode : %s, case : %s, bf-size : %s" % (topology, m, case, bf))

                    cd_lat = 0.0
                    link_usage = 0.0
                    if topology != '1337':

                        _data = read_data(data_dir[0], file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : case, 'mode' : m})

                        # bars  : link usage
                        # lines : avg. latency
                        cd_lat = float(_data['outcomes'].loc[_data['outcomes']['type'] == 'cd']['avg-latency'])

                        # SLM events are invariable for either mode
                        link_usage = _data['events']['slmi'].values + _data['events']['slmc'].values

                        # regarding MLM events:
                        #   - if mode is AML (m = 0), then we simply sum the avg. nr. of mlm events
                        #   - if mode is FB (m = 3), we must some the LLM events + the fwd events due 
                        #     to fallbacks, which are approx. the avg. latency of fallback events
                        if m == 0:
                            link_usage = (link_usage + _data['events']['mlm'].values) / 4.0
                        if m == 3:
                            __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']

                            print("link_usage = %s + %s + %s = %s" % 
                                (link_usage, _data['events']['llm'].values, 
                                __d['avg-latency'].values, 
                                link_usage + _data['events']['llm'].values + (__d['avg-latency'].values)))
                            link_usage = (link_usage + (__d['prob'].values * __d['avg-latency'].values)) / 4.0

                    else:
                        cd_lat = 4.0
                        link_usage = (4.0 / 4.0)

                    label = ['', '', '']
                    if setlabel == True:
                        label[0] = str(modes_str[str(bf)])
                        if m == 0:
                            label[1] = 'latency'

                    ax2.bar(
                        pos, 
                        link_usage,
                        color = colors[m], linewidth = 0.5, width = bar_width, label = label[0])

                    xx = 0.0
                    if (m == 0):
                        xx = bar_width

                    latency[m]['xx'].append(pos + xx)
                    latency[m]['yy'].append(cd_lat)

                    pos += bar_width

                    if (m == 0):
                        xbf[pos] = ('n' if (case[-1] == 'I') else (('n/2') if (case[-1] == 'H') else ('%s' % (case[-1]))))

                    if (m == 0) and (case == cases[metacase][0]):
                        xbf[pos + (1.25 * bar_width)] = '\n' + str(topology_labels[int(topology)])

                # only include lable once
                setlabel = False
                pos += (bar_width / 2.0)
            pos += (bar_width * 1.0)

            # lines hold avg. # of deliveries
            _label = {0 : '', 3 : ''}
            for _mode in latency:
                _bf = bf_sizes[_mode]

                if (t == 0):
                    _label[_mode] = modes_str[str(_bf)]

                ax1.plot(
                    latency[_mode]['xx'],
                    latency[_mode]['yy'],
                    linewidth = 0.0, color = 'black', linestyle = '-', markersize = 5, marker = markers[_bf], label = _label[_mode])
                # ax1.plot(
                #     latency[_mode]['xx'],
                #     latency[_mode]['yy'],
                #     linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = markers[_bf], label = _label[_mode])

        pos += (bar_width * 1.0)

    if scenario == 'cdn-lpm':

        # an horizontal shaded area, for fwd efficiency
        ax1.axhspan(4, 7, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.25)

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [4.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [5.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [6.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [7.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

    elif scenario == 'opportunistic':

        # an horizontal shaded area, for fwd efficiency
        ax1.axhspan(4, 5, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.25)

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [4.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [4.5] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [5.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

    # legends
    leg = []
    leg.append(ax1.legend(
        fontsize = 12, 
        ncol = 2, loc = 'upper left',
        handletextpad = 0.2, handlelength = 0.575, labelspacing = 0.2, columnspacing = 0.5, title = 'latency'))

    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 12)

    # axis labels
    ax1.set_xlabel("\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    # ax2.set_ylabel("avg. # of used links")
    ax1.set_ylabel("avg. latency")

    # avg used link axis
    ax2.set_xlim(0.0 - (bar_width / 2.0), pos - (2.0 * bar_width))
    ax2.set_ylim(0, 6)
    ax2.set_yticks([0, 1, 2])
    
    # fwd effic axis
    if scenario == 'cdn-lpm':
        ax1.set_ylim(-1, 11)
        ax1.set_yticks([4.0, 5.0, 6.0, 7.0])

    elif scenario == 'opportunistic':
        ax1.set_ylim(3, 6)
        ax1.set_yticks([4.0, 4.5, 5.0])

    # (1) avg. nr. of used links (bar chart)
    ax1 = fig.add_subplot(122)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(False)
    # (2) latency (line chart)
    ax2 = ax1.twinx()
    ax2.yaxis.grid(True, ls = 'dotted', lw = 0.50)
    ax1.set_zorder(ax2.get_zorder() + 1)
    ax1.patch.set_visible(False)

    pos = 0.0
    setlabel = True
    xbf = defaultdict()
    ax1.set_title("(b) w/ TPs", fontsize = 12)

    cases = OrderedDict()
    cstr = defaultdict()
    if scenario == 'opportunistic':

        metacases = ['{7,10}']
        _cases = {
            # '{7,10}' : ['S750221S1050221*172', 'S75022HS105022H*172', 'S75022IS105022I*172']
            '{7,10}' : ['S75022HS105022H*172', 'S75022IS105022I*172']
            }

    elif scenario == 'cdn-lpm':

        metacases = ['{4,5}']
        _cases = {
            '{2,3}' : ['S250221S350221*132', 'S25022HS35022H*132', 'S25022IS35022I*132'], 
            '{3,4}' : ['S350221S450221*142', 'S35022HS45022H*142', 'S35022IS45022I*142'], 
            # '{4,5}' : ['S450221S550221*152', 'S45022HS55022H*152', 'S45022IS55022I*152']
            '{4,5}' : ['S45022HS55022H*152', 'S45022IS55022I*152']
            }

    for metacase in metacases:
        cases[metacase] = _cases[metacase]

    latency = defaultdict()

    for t, topology in enumerate(topology_keys):

        latency[0] = defaultdict(list)
        latency[3] = defaultdict(list)

        for metacase in cases:
            for case in cases[metacase]:
                for m in modes:

                    bf = bf_sizes[m]
                    print("topology : %s, mode : %s, case : %s, bf-size : %s" % (topology, m, case, bf))

                    cd_lat = 0.0
                    link_usage = 0.0
                    if topology != '1337':

                        _data = read_data(data_dir[1], file_types = ['events', 'outcomes'], filters = {'topology' : topology, 'case' : case, 'mode' : m})

                        # bars  : link usage
                        # lines : avg. latency
                        cd_lat = float(_data['outcomes'].loc[_data['outcomes']['type'] == 'cd']['avg-latency'])

                        # SLM events are invariable for either mode
                        link_usage = _data['events']['slmi'].values + _data['events']['slmc'].values

                        # regarding MLM events:
                        #   - if mode is AML (m = 0), then we simply sum the avg. nr. of mlm events
                        #   - if mode is FB (m = 3), we must some the LLM events + the fwd events due 
                        #     to fallbacks, which are approx. the avg. latency of fallback events
                        if m == 0:
                            link_usage = (link_usage + _data['events']['mlm'].values) / 2.0
                        if m == 3:
                            __d = _data['outcomes'].loc[_data['outcomes']['type'] == 'hfd']

                            print("link_usage = %s + %s + %s = %s" % 
                                (link_usage, _data['events']['llm'].values, 
                                __d['avg-latency'].values, 
                                link_usage + _data['events']['llm'].values + (__d['avg-latency'].values)))
                            link_usage = (link_usage + (__d['prob'].values * __d['avg-latency'].values)) / 2.0

                    else:
                        cd_lat = 2.0
                        link_usage = (2.0 / 2.0)

                    label = ['', '', '']
                    if setlabel == True:
                        label[0] = str(modes_str[str(bf)])
                        if m == 0:
                            label[1] = 'latency'

                    ax2.bar(
                        pos, 
                        link_usage,
                        color = colors[m], linewidth = 0.5, width = bar_width, label = label[0])

                    xx = 0.0
                    if (m == 0):
                        xx = bar_width

                    latency[m]['xx'].append(pos + xx)
                    latency[m]['yy'].append(cd_lat)

                    pos += bar_width

                    if (m == 0):
                        xbf[pos] = ('n' if (case[-5] == 'I') else (('n/2') if (case[-5] == 'H') else ('%s' % (case[-5]))))

                    if (m == 0) and (case == cases[metacase][0]):
                        xbf[pos + (1.25 * bar_width)] = '\n' + str(topology_labels[int(topology)])

                # only include lable once
                setlabel = False
                pos += (bar_width / 2.0)
            pos += (bar_width * 1.0)

            # lines hold avg. # of deliveries
            _label = {0 : '', 3 : ''}
            for _mode in latency:
                _bf = bf_sizes[_mode]

                if (t == 0):
                    _label[_mode] = modes_str[str(_bf)]

                ax1.plot(
                    latency[_mode]['xx'],
                    latency[_mode]['yy'],
                    linewidth = 0.0, color = 'black', linestyle = '-', markersize = 5, marker = markers[_bf], label = _label[_mode])
                # ax1.plot(
                #     latency[_mode]['xx'],
                #     latency[_mode]['yy'],
                #     linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = markers[_bf], label = _label[_mode])

        pos += (bar_width * 1.0)

    if scenario == 'cdn-lpm':
        # an horizontal shaded area, for fwd efficiency
        ax1.axhspan(2, 4, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.25)

        # ax1.plot(
        #     [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        #     [1.0] * 2,
        #     linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [2.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [3.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [4.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

    elif scenario == 'opportunistic':

        # an horizontal shaded area, for fwd efficiency
        ax1.axhspan(2, 2.5, linewidth = 0.0, facecolor = '#bebebe', alpha = 0.25)

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [2.0] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        ax1.plot(
            [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
            [2.5] * 2,
            linewidth = 0.5, color = 'black', linestyle = ':')

        # ax1.plot(
        #     [0.0 - (bar_width / 2.0), pos - (1.0 * bar_width)],
        #     [3.0] * 2,
        #     linewidth = 0.5, color = 'black', linestyle = ':')


    # legends
    leg = []
    leg.append(ax2.legend(
        fontsize = 12, 
        ncol = 1, loc = 'upper right', title = 'link usage exc.',
        handletextpad = 0.2, handlelength = 1.0, labelspacing = 0.2, columnspacing = 0.5))
    # set legend title fontsize to 10
    for l in leg:
        plt.setp(l.get_title(), fontsize = 12)

    # axis labels
    ax1.set_xlabel("\n")
    ax1.set_xticks(xbf.keys())
    ax1.set_xticklabels(xbf.values())

    ax2.set_ylabel("link usage excess")
    # ax1.set_ylabel("avg. latency")

    fig.text(0.50, 0.05, "# of caches around request source\ntopology", ha = 'center')

    # avg used link axis
    ax2.set_xlim(0.0 - (bar_width / 2.0), pos - (2.0 * bar_width))
    ax2.set_ylim(0, 6)
    ax2.set_yticks([0, 1, 2])
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 6, 1)])
    # ax2.set_ylim(1, 10000)
    # ax2.set_yscale("log", nonposy = 'clip')
    # ax2.set_yticks([10**(x) for x in np.arange(0, 5, 1)])
    # ax2.tick_params(axis = 'y', which = 'minor', bottom = 'off')
    
    # fwd effic axis
    if scenario == 'cdn-lpm':
        ax1.set_ylim(-2, 8)
        ax1.set_yticks([2.0, 3.0, 4.0])
    elif scenario == 'opportunistic':
        ax1.set_ylim(1, 4)
        ax1.set_yticks([2.0, 2.5])

    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, ("%s-usage.pdf" % (scenario))), bbox_inches='tight', format = 'pdf')
