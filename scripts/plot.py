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

matplotlib.rcParams.update({'font.size': 16})

FONTSIZE_LEGEND = 14

def extract_data(data_dir):

    data = defaultdict(OrderedDict)

    for file_name in sorted(glob.glob(os.path.join(data_dir, '*.tsv'))):

        file_type = file_name.split(".")[0].split("/")[-1]
        file_label = file_name.split(".")[1]

        print("type = %s, label = %s" % (file_type, file_label))

        data[file_type][file_label] = pd.read_csv(file_name, sep = "\t")
        data[file_type][file_label] = data[file_type][file_label].convert_objects(convert_numeric = True)

    return data

def plot_base(data, subcase):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    event_avgs = defaultdict(list)
    path_values = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Single link', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe', 'white']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black']
    path_styles = ['-', '--']
    path_markers = ['o', 'v']

    # keeps the x-axis labels
    entry_lengths = []

    for key in data:

        for sub_key in data[key]:

            if key == "events":

                # strip any 'F0' prefix
                entry_lengths.append(sub_key.upper().lstrip("F").lstrip("0"))
                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(key)
                print(sub_key)
                print(event_probs)

                # we want to present the bars is a diff. order than the 
                # event codes, hence the use of a dict
                for k, i in {0:0, 1:3, 2:1, 3:2}.iteritems():
                    v = (event_probs[i])
                    event_avgs[k].append(v)

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()
                path_slots = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].count()

                print(key)
                print(sub_key)
                print(path_probs)
                print(path_slots)

                for k, i in {0:0, 1:1}.iteritems():

                    if subcase == "flood":
                        # FIXME: there are only 6 possible slots for wrong 
                        # deliveries, path_slots[1] is equal to 15 though
                        if i == 1:
                            path_slots[i] = 6.0

                        v = (path_probs[i] / path_slots[i])

                    else:

                        v = path_probs[i]

                    path_values[k].append(v)

    # always add the 'ideal' case: this represents the avg. number of events 
    # for a completely correct forwarding case:
    #   -# single egress iface: 6 hops
    #   -# local iface: 1 hop (the correct content source)
    #   -# no packet drops or multiple matches
    #   -# 1.0 correct deliv. prob
    #   -# 0.0 wrong deliv. prob
    for k, v in {0: 0.0, 1: 6.0, 2: 0.0, 3: 1.0, 4: 1.0, 5: 0.0}.iteritems():

        if k < 4:
            event_avgs[k].append(v)
        else:
            path_values[k - 4].append(v)

    entry_lengths.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)

    m = -(float(len(event_labels)) / 2.0)
    bar_width = 0.30
    
    for l in np.arange(len(event_labels)):
        ax1.bar(np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width), np.array(event_avgs[l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = event_labels[l])
        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.01, 1000.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|R|=15, 192 bit BF)")
    ax1.set_xlabel("Fwd. entry size")
    ax1.set_ylabel("Avg. nr. of events")
    ax1.set_xticks(np.arange(1, (2 * len(entry_lengths)), step = 2), entry_lengths)
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(np.arange(1, (2 * len(entry_lengths)), step = 2), np.array(path_values[l]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], marker = path_markers[l], label = path_labels[l])

    ax2.set_ylim(0.0, 2.5)
    ax2.set_ylabel("Path outcome prob.")
    ax2.set_xticks(np.arange(1, (2 * len(entry_lengths)), step = 2))
    ax2.set_xticklabels(entry_lengths)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("base-" + subcase + ".pdf", bbox_inches='tight', format = 'pdf')

def plot_bf_sizes(data):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    event_avgs = defaultdict(list)
    path_values = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Single link', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe', 'white']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black']
    path_styles = ['-', '--']
    path_markers = ['o', 'v']

    # keeps the x-axis labels
    bf_sizes = []

    for key in data:

        for sub_key in data[key]:

            if key == "events":

                # strip any 'F0' prefix
                bf_sizes.append(sub_key.upper().lstrip("BF").lstrip("0"))
                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(key)
                print(sub_key)
                print(event_probs)

                # we want to present the bars is a diff. order than the 
                # event codes, hence the use of a dict
                for k, i in {0:0, 1:3, 2:1, 3:2}.iteritems():
                    v = (event_probs[i])
                    event_avgs[k].append(v)

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()
                path_slots = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].count()

                for k, i in {0:0, 1:1}.iteritems():

                    # FIXME: there are only 6 possible slots for wrong 
                    # deliveries, path_slots[1] is equal to 15 though
                    if i == 1:
                        path_slots[i] = 6.0

                    v = (path_probs[i] / path_slots[i])
                    path_values[k].append(v)

    # always add the 'ideal' case: this represents the avg. number of events 
    # for a completely correct forwarding case:
    #   -# single egress iface: 6 hops
    #   -# local iface: 1 hop (the correct content source)
    #   -# no packet drops or multiple matches
    #   -# 1.0 correct deliv. prob
    #   -# 0.0 wrong deliv. prob
    for k, v in {0: 0.0, 1: 6.0, 2: 0.0, 3: 1.0, 4: 1.0, 5: 0.0}.iteritems():

        if k < 4:
            event_avgs[k].append(v)
        else:
            path_values[k - 4].append(v)

    bf_sizes.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)

    m = -(float(len(event_labels)) / 2.0)
    bar_width = 0.30
    
    for l in np.arange(len(event_labels)):
        ax1.bar(np.arange(1, (2 * len(bf_sizes)), step = 2) + (m * bar_width), np.array(event_avgs[l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = event_labels[l])
        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.01, 1000.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|R|=15, |F|=1)")
    ax1.set_xlabel("BF sizes (in bit)")
    ax1.set_ylabel("Avg. nr. of events")
    ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(np.arange(1, (2 * len(bf_sizes)), step = 2), np.array(path_values[l]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], marker = path_markers[l], label = path_labels[l])

    ax2.set_ylim(0.0, 2.5)
    ax2.set_ylabel("Path outcome prob.")
    ax2.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2))
    ax2.set_xticklabels(bf_sizes)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("bf-sizes.pdf", bbox_inches='tight', format = 'pdf')

def plot_req_sizes(data):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    event_avgs = defaultdict(list)
    path_values = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Single link', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe', 'white']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black']
    path_styles = ['-', '--']
    path_markers = ['o', 'v']

    # keeps the x-axis labels
    req_sizes = []

    for key in data:

        for sub_key in data[key]:

            if key == "events":

                # strip any 'F0' prefix
                req_sizes.append(sub_key.upper().lstrip("R").lstrip("0"))
                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(key)
                print(sub_key)
                print(event_probs)

                # we want to present the bars is a diff. order than the 
                # event codes, hence the use of a dict
                for k, i in {0:0, 1:3, 2:1, 3:2}.iteritems():
                    v = (event_probs[i])
                    event_avgs[k].append(v)

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()
                path_slots = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].count()

                for k, i in {0:0, 1:1}.iteritems():

                    # FIXME: there are only 6 possible slots for wrong 
                    # deliveries, path_slots[1] is equal to 15 though
                    if i == 1:
                        path_slots[i] = 6.0

                    v = (path_probs[i] / path_slots[i])
                    path_values[k].append(v)

    # always add the 'ideal' case: this represents the avg. number of events 
    # for a completely correct forwarding case:
    #   -# single egress iface: 6 hops
    #   -# local iface: 1 hop (the correct content source)
    #   -# no packet drops or multiple matches
    #   -# 1.0 correct deliv. prob
    #   -# 0.0 wrong deliv. prob
    for k, v in {0: 0.0, 1: 6.0, 2: 0.0, 3: 1.0, 4: 1.0, 5: 0.0}.iteritems():

        if k < 4:
            event_avgs[k].append(v)
        else:
            path_values[k - 4].append(v)

    req_sizes.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5,4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)

    m = -(float(len(event_labels)) / 2.0)
    bar_width = 0.30
    
    for l in np.arange(len(event_labels)):
        ax1.bar(np.arange(1, (2 * len(req_sizes)), step = 2) + (m * bar_width), np.array(event_avgs[l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = event_labels[l])
        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.01, 1000.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("Req. size")
    ax1.set_ylabel("Avg. nr. of events")
    ax1.set_xticks(np.arange(1, (2 * len(req_sizes)), step = 2), req_sizes)
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(np.arange(1, (2 * len(req_sizes)), step = 2), np.array(path_values[l]), linewidth = 1.5, linestyle = path_styles[l], color = path_colors[l], marker = path_markers[l], label = path_labels[l])

    ax2.set_ylim(0.0, 2.5)
    ax2.set_ylabel("Path outcome prob.")
    ax2.set_xticks(np.arange(1, (2 * len(req_sizes)), step = 2))
    ax2.set_xticklabels(req_sizes)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("req-sizes.pdf", bbox_inches='tight', format = 'pdf')

def plot_cdn_latencies(data):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    latency_avgs    = defaultdict(list)
    event_avgs      = defaultdict(defaultdict)
    path_values     = defaultdict(defaultdict)

    # list of exclusion words
    exclusions = ['192']

    # labels and colors for bars (events)
    event_labels = ['Correct deliv.', 'Wrong deliv.']
    event_colors = ['lightgrey', '#708090']
    event_hatches = ['x', 'o', '*']
    # labels and styles for lines (latency probs)
    latency_labels = ['Flood', 'Random', 'Fallback', 'Ideal']
    latency_colors = ['black', 'black', 'black', 'black']
    latency_styles = ['-', '--', '-.', ':']
    latency_markrs = ['x', 'v', '^', 'o']
    latency_markrs_size = [14, 10, 10, 8]

    bf_sizes = []

    for key in data:

        for sub_key in data[key]:

            if "192" in sub_key:
                continue

            # filenames follow the format <type>.BF<bf-size>-<mode>.tsv
            # extract the mode
            mode = int(sub_key.split("-", 1)[1])

            if key == "events":

                # strip any 'BF0' prefix
                bf_size = sub_key.split("-", 1)[0].upper().lstrip("BF").lstrip("0")

                if bf_size not in bf_sizes:
                    bf_sizes.append(bf_size)

                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                if mode == 2:
                    print(sub_key)
                    print(event_probs)

                for i in [0, 1]:

                    if i not in event_avgs[mode]:
                        event_avgs[mode][i] = []

                    # if 'Fallback', then we will have deliveries from 2 causes:
                    #   -# fallbacks: when MLM event happens
                    #   -# normal RID forwarding behavior: when LLM events happen
                    #
                    # the sum of the 2 should ideally be 1, but let's present 
                    # the value of the sum, just in case
                    if mode == 2:

                        if i == 0:
                            v = event_probs[1] + event_probs[2]
                        else:
                            v = event_probs[2]

                        event_avgs[mode][i].append(v)

                    else:
                        # otherwise, only the LLM event (local link matches) 
                        # stats matter
                        event_avgs[mode][i].append(event_probs[2])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                if mode == 2:
                    print(sub_key)
                    print(path_probs)

                for i in [0, 1]:

                    # if mode is 'Fallback', deliveries are guaranteed to be 
                    # correct
                    if mode == 2:

                        if i == 0:
                            v = 1.0
                        else:
                            v = 0.0

                    else:
                        v = path_probs[i]

                    if i not in path_values[mode]:
                        path_values[mode][i] = []

                    path_values[mode][i].append(v)

                avg_latency = 0.0
                latency_probs = data[key][sub_key].groupby(by = ["LATENCY"])["PROB"].sum()

                # avg. latency is generally calculated as sum(latency * pron), 
                # but if the MMH_MODE is 'Flood', we simply take the min() 
                # of the keys for which prob > 0.0
                if mode == 0:

                    avg_latency = float(min(latency_probs.loc[latency_probs > 0.0].index.tolist()))

                else:

                    for lat, prob in latency_probs.iteritems():
                        avg_latency += (lat * prob)

                latency_avgs[mode].append(avg_latency)

    # scale the LLM avg. values by correctness
    for mode in event_avgs:
        for status in event_avgs[mode]:
            for i in np.arange(len(event_avgs[mode][status])):

                if mode == 0:
                    event_avgs[mode][status][i] = path_values[mode][status][i]
                else:
                    event_avgs[mode][status][i] = event_avgs[mode][status][i] * path_values[mode][status][i]

    # always add the 'ideal' case for comparison:
    #   -# LLM: 1
    #   -# avg. latency: 3
    # event_avgs[0][0].append(0.0)
    # event_avgs[0][1].append(0.0)
    # event_avgs[1][0].append(1.0)
    # event_avgs[1][1].append(0.0)
    # event_avgs[2][0].append(0.0)
    # event_avgs[2][1].append(0.0)

    # for m in np.arange(3):
    #     latency_avgs[m].append(3.0)

    for l in latency_avgs[0]:
        latency_avgs[3].append(l)

    # bf_sizes.append('Ideal')

    print("latency")
    print(latency_avgs)

    print("events")
    for k in event_avgs:
        print(k)
        print(event_avgs[k])

    print(bf_sizes)

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    
    # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 3

    # assumes the inter bar group space is half a bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    show_legend = True
    xticks = []
    for mode in np.arange(3):

        if mode > 0:
            show_legend = False

        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(bf_sizes)), step = 2) + (m * bar_width), np.array(event_avgs[mode][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.0001, 1000.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("BF sizes (in bit)")
    ax1.set_ylabel("Avg. nr. of deliv.")
    ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.set_yticks([0.0001, 0.001, 0.01, 0.1, 1.0, 10.0])
#    ax1.set_yticklabels(['$10^{-4}$', '10^-3', '10^-2', '10^-1', '1', '10', '', ''])

    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(latency_labels)):
        ax2.plot(np.arange(1, (2 * len(bf_sizes)), step = 2), np.array(latency_avgs[l]), linewidth = 1.5, color = latency_colors[l], linestyle = latency_styles[l], marker = latency_markrs[l], markersize = latency_markrs_size[l], label = latency_labels[l])

    ax2.set_ylim(0.0, 14.0)
    ax2.set_ylabel("Avg. latency")
#    ax2.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2))
#    xtick_labels = bf_sizes
    xticks = [0.4, 1.0, 1.6, 2.4, 3.0, 3.6, 4.4, 5.0, 5.6]
    # xtick_labels = ['FL', 'R\n256', 'FB', 'FL', 'R\n512', 'FB', 'FL', 'R\n1024', 'FB']
    xtick_labels = ['FL', 'R\n256', 'FB', '', '\n512', '', '', '\n1024', '']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
    # ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("cache-cdn-latencies.pdf", bbox_inches='tight', format = 'pdf')

def plot_opportunistic_latencies(data):

    latency_avgs    = defaultdict(list)
    event_avgs      = defaultdict(defaultdict)
    path_values     = defaultdict(defaultdict)

    # labels and colors for bars (events)
    event_labels = ['Correct deliv.', 'Wrong deliv.']
    event_colors = ['lightgrey', '#708090']
    event_hatches = ['x', 'o', '*']
    # labels and styles for lines (latency probs)
    latency_labels = ['|F|:7', '|F|:15', 'Ideal']
    latency_colors = ['black', 'black', 'black']
    latency_styles = ['-', '--', '-.']
    latency_markrs = ['x', 'v', 'o']
    latency_markrs_size = [14, 10, 8]

    cache_distances = []
    bf_sizes = []
    probs_sum = defaultdict(list)

    for key in data:
        for sub_key in data[key]:

            # filenames follow the format: 
            # <type>.<cache-distance>-<bf-size>:L<local-entry-size>.tsv

            # extract local entry size
            entry_size = int(sub_key.split(":", 1)[1].upper().lstrip("L").lstrip("0"))

            if key == "events":

                # extract the cache dist and bf size
                cache_distance = sub_key.split(":", 1)[0]
                # bf_size = int((sub_key.split(":", 1)[0]).split("-", 1)[1].upper().lstrip("BF").lstrip("0"))

                if cache_distance not in cache_distances:
                    cache_distances.append(cache_distance)

                # if bf_size not in bf_sizes:
                #     bf_sizes.append(bf_size)

                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(sub_key)
                print(event_probs)

                for i, e in {0:2, 1:2}.iteritems():

                    if i not in event_avgs[entry_size]:
                        event_avgs[entry_size][i] = []
                    
                    event_avgs[entry_size][i].append(event_probs[e])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                print(sub_key)
                print(path_probs)

                p = 0.0
                for i in [0, 1]:

                    v = path_probs[i]

                    if i not in path_values[entry_size]:
                        path_values[entry_size][i] = []

                    path_values[entry_size][i].append(v)
                    p += v

                probs_sum[entry_size].append(p)

                avg_latency = 0.0
                latency_probs = data[key][sub_key].groupby(by = ["LATENCY"])["PROB"].sum()

                for lat, prob in latency_probs.iteritems():
                    avg_latency += (lat * prob)

                latency_avgs[entry_size].append(avg_latency)

    print("probs")
    for k in path_values:
        for i in [0, 1]:
            print(k)
            print(path_values[k][i])
            print(probs_sum[k])
            path_values[k][i] = [ a / b for a, b in zip(path_values[k][i], probs_sum[k]) ]
            print(path_values[k][i])
            
    # scale the LLM avg. values by correctness
    for entry_size in event_avgs:
        for status in event_avgs[entry_size]:
            for i in np.arange(len(event_avgs[entry_size][status])):

                v = event_avgs[entry_size][status][i] * path_values[entry_size][status][i]

                # ok, this is cheating a bit: the idea is to show a bit 
                # of a bar (representing 0.0), just for the reader to know that this 
                # isn't an error
                if v < 0.005:
                    v = 0.005

                event_avgs[entry_size][status][i] = v

    # ideal latencies (added manually)
    latency_avgs[16].append(2.0)
    latency_avgs[16].append(3.0)
    latency_avgs[16].append(3.0)
    latency_avgs[16].append(3.0)

    print("latency")
    print(latency_avgs)

    print("events")
    for k in event_avgs:
        print(k)
        print(event_avgs[k])

    print("cache distances")
    print(cache_distances)

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    
    # we have 2 local entry sizes and 2 correctness values per size. hence 2 
    # groups of bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 2

    # assumes the inter bar group space is 1 bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_min = 0.0

    bar_group_offset = np.arange(1, (2 * len(cache_distances)), step = 2)[0] + (m * bar_width)
    bar_group_width = (bar_group_num * bar_group_size + 1) * bar_width

    x_max = np.arange(1, (2 * len(cache_distances)), step = 2)[0] + (m * bar_width)
    x_max += np.arange(1, (2 * len(cache_distances)), step = 2)[-1] + ((bar_group_num * bar_group_size + 1) / 2.0) * bar_width

    print(bar_group_size)
    print(bar_group_num)
    print(bar_group_offset)
    print(bar_group_width)
    print("[%d, %d]" % (x_min, x_max))

    show_legend = True
    for entry_size in [7, 15]:

        # only show the legend for 1 entry size
        if entry_size > 7:
            show_legend = False

        # print the bars for each |L|, for both cache distances
        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(cache_distances)), step = 2) + (m * bar_width), np.array(event_avgs[entry_size][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        m += 1.0

#    ax1.set_yscale('log')
    ax1.set_ylim(0.0, 1.75)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("Cache distances")
    ax1.set_ylabel("Avg. nr. of deliv.")
    # ax1.set_xticks(np.arange(1, (2 * len(cache_distances)), step = 2), cache_distances)
    ax1.set_yticks(np.arange(0.0, 1.75, step = 0.25))
    ax1.set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0', '1.25', ''])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for k, l in {0:7, 1:15, 2:16}.iteritems():

        ax2.plot(np.arange(1, (2 * len(cache_distances)), step = 2), np.array(latency_avgs[l]), linewidth = 1.5, color = latency_colors[k], linestyle = latency_styles[k], marker = latency_markrs[k], markersize = latency_markrs_size[k], label = latency_labels[k])

    ax2.set_ylim(0.0, 14.0)
    ax2.set_ylabel("Avg. latency")

    xticks = [0.2, 0.7, 1.0, 1.3, 2.7, 3.0, 3.3, 4.7, 5.0, 5.3, 6.7, 7.0, 7.3,]
    xtick_labels = ['|F|:', '7', '\n2', '15', '7', '\n3', '15', '7', '\n3$\,$(FB)', '15', '7', '\n3$\,$(512)', '15']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    # ax2.set_xlim(x_min, x_max)

    ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
    # ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("cache-opportunistic-latencies.pdf", bbox_inches='tight', format = 'pdf')

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
         help = """the case you want to output. e.g. 'base', 'bf-sizes', etc.""")
    parser.add_argument(
        "--subcase", 
         help = """the sub-case you want to output. e.g. 'flood' or 'random' for 'base'.""")

    args = parser.parse_args()

    # quit if a dir w/ causality files hasn't been provided
    if not args.data_dir:
        sys.stderr.write("""%s: [ERROR] please supply a data dir!\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)

    # if an output dir is not specified, use data-dir
    if not args.output_dir:
        args.output_dir = args.data_dir

    # extract the data from all files in the data dir
    data = extract_data(args.data_dir)

    if args.case == 'base':

        if args.subcase in ['flood', 'random']:
            plot_base(data, args.subcase)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase ('flood' or 'random').\n""" % sys.argv[0]) 
            parser.print_help()
            sys.exit(1)            

    elif args.case == 'bf-sizes':
        plot_bf_sizes(data)
    elif args.case == 'req-sizes':
        plot_req_sizes(data)
    elif args.case == 'cache-pro-active':
        plot_cdn_latencies(data)
    elif args.case == 'cache-opportunistic':
        plot_opportunistic_latencies(data)
    else:
        sys.stderr.write("""%s: [ERROR] please supply a valid case\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)