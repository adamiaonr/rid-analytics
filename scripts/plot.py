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

def interleave(a, b):
    c = list(zip(a, b))
    return [elt for sublist in c for elt in sublist]

def extract_data(data_dir):

    data = defaultdict(OrderedDict)

    for file_name in sorted(glob.glob(os.path.join(data_dir, '*.tsv'))):

        file_type = file_name.split(".")[0].split("/")[-1]
        file_label = file_name.split(".")[1]

        print("filename = %s, type = %s, label = %s" % (file_name, file_type, file_label))

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

def plot_global_topologies_entry_sizes(data):

    event_avgs = defaultdict(defaultdict)
    path_values = defaultdict(defaultdict)
    delivery_probs = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black', 'black', 'black']
    path_styles = ['-', '--', '-.', ':']
    path_markrs = ['o', 'v', '^', 'o']
    path_markrs_colors = ['white', 'white', 'black', 'black']
    path_markrs_size = [10, 10, 10, 8]

    # x-axis labels
    entry_lengths = []

    for key in data:

        for sub_key in data[key]:

            # first letter of sub_key is topology type
            topology_type = sub_key.split("-", 1)[0][0]

            if key == "events":

                entry_length = int(sub_key.split("-", 1)[0].upper().lstrip("H").lstrip("D").lstrip("0"))
                if entry_length not in entry_lengths:
                    entry_lengths.append(entry_length)

                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:2}.iteritems():

                    if k not in event_avgs[topology_type]:
                        event_avgs[topology_type][k] = []

                    event_avgs[topology_type][k].append(event_probs[i])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:4}.iteritems():

                    if k not in path_values[topology_type]:
                        path_values[topology_type][k] = []

                    if i not in path_probs:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[topology_type][k].append(v)

    for i in np.arange(len(path_values["D"][0])):
        for s in [0, 1]:
            for topology_type in ['H', 'D']:
                delivery_probs[s].append(path_values[topology_type][s][i])

    # FIXME: remove the ideal probabilities
    # for k, v in {0: 1.0, 1: 0.0}.iteritems():
    #     delivery_probs[k].append(v)
    #     delivery_probs[k].append(v)

    print("DELIVERY PROBS")
    print(delivery_probs)

    # # always add the 'ideal' case: this represents the avg. number of events 
    # # for a completely correct forwarding case:
    # #   -# single egress iface: 6 hops
    # #   -# local iface: 1 hop (the correct content source)
    # #   -# no packet drops or multiple matches
    # #   -# 1.0 correct deliv. prob
    # #   -# 0.0 wrong deliv. prob
    # for k, v in {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 0.0}.iteritems():
    #     for mode in [0, 1]:
    #         if k < 3:
    #             event_avgs[mode][k].append(v)
    #         else:
    #             path_values[mode][k - 3].append(v)

    # entry_lengths.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.yaxis.grid(True)
    
    # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 2

    # assumes the inter bar group space is half a bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_pos = defaultdict(list)

    show_legend = True
    for topology_type in ["H", "D"]:

        if topology_type == "D":
            show_legend = False

        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width), np.array(event_avgs[topology_type][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        # positions of centers of bar groups
        x_pos[topology_type] = np.arange(1, (2 * len(entry_lengths)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    print("SHIIIIAAAT")
    print(interleave(x_pos["H"], x_pos["D"]))
    print(delivery_probs)

    ax1.set_yscale('log')
    ax1.set_ylim(0.0001, 1000000.0)

    ax1.set_xlabel("Fwd. entry size")
    ax1.set_ylabel("Avg. nr. of events.")
    # ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.set_yticks([0.0001, 0.01, 1.0, 100.0])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(interleave(x_pos["H"], x_pos["D"]), np.array(delivery_probs[l]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], markersize = path_markrs_size[l], marker = path_markrs[l], label = path_labels[l])
        ax2.plot(interleave(x_pos["H"], x_pos["D"]), np.array(delivery_probs[l]), linewidth = 1.5, color = path_markrs_colors[l], linestyle = '', markersize = path_markrs_size[l] - 3, marker = path_markrs[l], label = '')

    ax2.set_ylim(-0.25, 10.00)
    ax2.set_ylabel("Avg. nr. of deliveries")

    xticks = [0.6, 1.0, 1.4, 2.6, 3.0, 3.4, 4.6, 5.0, 5.4, 6.6, 7.0, 7.4]
    xtick_labels = ['H', '\n1', 'NH', '', '\n5', '', '', '\n10', '', '', '\n15', '']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("global-topologies-entry-size.pdf", bbox_inches='tight', format = 'pdf')

def plot_global_topologies_bf_sizes(data):

    event_avgs = defaultdict(defaultdict)
    path_values = defaultdict(defaultdict)
    delivery_probs = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black', 'black', 'black']
    path_styles = ['-', '--', '-.', ':']
    path_markrs = ['o', 'v', '^', 'o']
    path_markrs_colors = ['white', 'white', 'black', 'black']
    path_markrs_size = [10, 10, 10, 8]

    # x-axis labels
    bf_sizes = []

    for key in data:

        for sub_key in data[key]:

            topology_type = sub_key.split("-", 1)[1]

            if key == "events":

                bf_size = int(sub_key.split("-", 1)[0].upper().lstrip("BF").lstrip("0"))
                if bf_size not in bf_sizes:
                    bf_sizes.append(bf_size)

                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:2}.iteritems():

                    if k not in event_avgs[topology_type]:
                        event_avgs[topology_type][k] = []

                    event_avgs[topology_type][k].append(event_probs[i])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:4}.iteritems():

                    if k not in path_values[topology_type]:
                        path_values[topology_type][k] = []

                    if i not in path_probs:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[topology_type][k].append(v)

    # for i in np.arange(len(path_values[1][0])):
    #     path_values[1][0][i] = path_values[1][0][i] - path_values[1][1][i] - path_values[1][2][i]

    # print("DA THINGA")
    # print(path_values[1][2])
    # print(path_values[1][1])
    # print(path_values[1][0])

    # for i in np.arange(len(event_avgs[1][2])):
    #     event_avgs[1][2][i] = event_avgs[1][2][i] - path_values[1][1][i] - path_values[1][2][i]

    for i in np.arange(len(path_values["D"][0])):
        for s in [0, 1]:
            for topology_type in ["H", "D"]:
                delivery_probs[s].append(path_values[topology_type][s][i])

    # for k,v in {0: 1.0, 1: 0.0}.iteritems():
    #     delivery_probs[k].append(v)
    #     delivery_probs[k].append(v)

    # print("DELIVERY PROBS")
    # print(delivery_probs)

    # # always add the 'ideal' case: this represents the avg. number of events 
    # # for a completely correct forwarding case:
    # #   -# single egress iface: 6 hops
    # #   -# local iface: 1 hop (the correct content source)
    # #   -# no packet drops or multiple matches
    # #   -# 1.0 correct deliv. prob
    # #   -# 0.0 wrong deliv. prob
    # for k, v in {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 0.0}.iteritems():
    #     for mode in [0, 1]:
    #         if k < 3:
    #             event_avgs[mode][k].append(v)
    #         else:
    #             path_values[mode][k - 3].append(v)

    # bf_sizes.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.yaxis.grid(True)
    
    # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 2

    # assumes the inter bar group space is half a bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_pos = defaultdict(list)

    show_legend = True
    for topology_type in ["H", "D"]:

        if topology_type == "D":
            show_legend = False

        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(bf_sizes)), step = 2) + (m * bar_width), np.array(event_avgs[topology_type][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        x_pos[topology_type] = np.arange(1, (2 * len(bf_sizes)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    print(interleave(x_pos["H"], x_pos["D"]))

    ax1.set_yscale('log')
    ax1.set_ylim(0.000001, 100000.0)

    ax1.set_xlabel("BF sizes (bit)")
    ax1.set_ylabel("Avg. nr. of events.")
    # ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.set_yticks([0.000001, 0.0001, 0.01, 1.0, 10.0])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(interleave(x_pos["H"], x_pos["D"]), np.array(delivery_probs[l]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], markersize = path_markrs_size[l], marker = path_markrs[l], label = path_labels[l])
        ax2.plot(interleave(x_pos["H"], x_pos["D"]), np.array(delivery_probs[l]), linewidth = 1.5, color = path_markrs_colors[l], linestyle = '', markersize = path_markrs_size[l] - 3, marker = path_markrs[l], label = '')

    ax2.set_ylim(-0.05, 2.00)
    ax2.set_ylabel("Avg. nr. of deliveries")

    xticks = [0.6, 1.0, 1.4, 2.6, 3.0, 3.4, 4.6, 5.0, 5.4]
    xtick_labels = ['H', '\n256', 'NH', '', '\n512', '', '', '\n1024', '']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    # ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("global-topologies-bf-size.pdf", bbox_inches='tight', format = 'pdf')

def plot_prefix_exclusion(data):

    event_avgs = defaultdict(defaultdict)
    path_values = defaultdict(defaultdict)
    delivery_probs = defaultdict(defaultdict)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe']
    # labels and styles for lines (path probs)
    path_labels = ['Cache pos.: 3.1, Annc. Radius: 2', 'CP: 3.1, AR: 3', 'CP: 3.3, AR: 3', 'CP: 3.3, AR: 4']
    path_colors = ['black', 'black', 'black', 'black']
    path_styles = ['-', '-', '--', '--']
    path_markrs = ['o', 'v', '^', 'o']
    path_markrs_colors = ['white', 'white', 'white', 'white']
    path_markrs_size = [10, 10, 10, 8]

    # x-axis labels
    req_sizes = []

    for key in data:

        for sub_key in data[key]:

            if sub_key.split("-", 1)[1] == "0":
                continue

            cache_distance = sub_key.split("-", 1)[0].split(":")[1]
            radius = sub_key.split("-", 1)[0].split(":")[0]

            if key == "events":

                req_size = int(sub_key.split("-", 1)[0].split(":")[2].upper().lstrip("0"))
                if req_size not in req_sizes:
                    req_sizes.append(req_size)

                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:2}.iteritems():

                    if radius not in event_avgs[cache_distance]:
                        event_avgs[cache_distance][radius] = dict()

                    if k not in event_avgs[cache_distance][radius]:
                        event_avgs[cache_distance][radius][k] = []

                    event_avgs[cache_distance][radius][k].append(event_probs[i])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:4}.iteritems():

                    if radius not in path_values[cache_distance]:
                        path_values[cache_distance][radius] = dict()

                    if k not in path_values[cache_distance][radius]:
                        path_values[cache_distance][radius][k] = []

                    if i not in path_probs:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[cache_distance][radius][k].append(v)

    # for i in np.arange(len(path_values[1][0])):
    #     path_values[1][0][i] = path_values[1][0][i] - path_values[1][1][i] - path_values[1][2][i]

    # print("DA THINGA")
    # print(path_values[1][2])
    # print(path_values[1][1])
    # print(path_values[1][0])

    # for i in np.arange(len(event_avgs[1][2])):
    #     event_avgs[1][2][i] = event_avgs[1][2][i] - path_values[1][1][i] - path_values[1][2][i]

    for i in np.arange(len(path_values["3"]["SF"][0])):
        for cache_distance in ["1", "3"]:
            for radius in ["S", "SF"]:
                for s in [0, 1]:

                    if radius not in delivery_probs[cache_distance]:
                        delivery_probs[cache_distance][radius] = dict()

                    if s not in delivery_probs[cache_distance][radius]:
                        delivery_probs[cache_distance][radius][s] = []

                    delivery_probs[cache_distance][radius][s].append(path_values[cache_distance][radius][s][i])

    # for k,v in {0: 1.0, 1: 0.0}.iteritems():
    #     delivery_probs[k].append(v)
    #     delivery_probs[k].append(v)

    # print("DELIVERY PROBS")
    # print(delivery_probs)

    # # always add the 'ideal' case: this represents the avg. number of events 
    # # for a completely correct forwarding case:
    # #   -# single egress iface: 6 hops
    # #   -# local iface: 1 hop (the correct content source)
    # #   -# no packet drops or multiple matches
    # #   -# 1.0 correct deliv. prob
    # #   -# 0.0 wrong deliv. prob
    # for k, v in {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 0.0}.iteritems():
    #     for mode in [0, 1]:
    #         if k < 3:
    #             event_avgs[mode][k].append(v)
    #         else:
    #             path_values[mode][k - 3].append(v)

    # bf_sizes.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.yaxis.grid(True)
    
    # # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # # bars, each group with 2 bars.
    # bar_group_size = len(event_labels)
    # bar_group_num = 2

    # # assumes the inter bar group space is half a bar. also, for n bar groups 
    # # we have n - 1 inter bar group spaces
    # m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    # bar_width = 0.20

    # x_pos = defaultdict(list)

    # show_legend = True
    # for topology_type in ["H", "D"]:

    #     if topology_type == "D":
    #         show_legend = False

    #     for l in np.arange(len(event_labels)):

    #         if show_legend:
    #             leg = event_labels[l]
    #         else:
    #             leg = None

    #         ax1.bar(np.arange(1, (2 * len(bf_sizes)), step = 2) + (m * bar_width), np.array(event_avgs[topology_type][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
    #         m += 1.0

    #     x_pos[topology_type] = np.arange(1, (2 * len(bf_sizes)), step = 2) + ((m - 1.5) * bar_width)
    #     m += 1.0

    # print(interleave(x_pos["H"], x_pos["D"]))

    # ax1.set_yscale('log')
    # ax1.set_ylim(0.000001, 100000.0)

    ax1.set_xlabel("Nr. of prefix lengths ($|R|_{max}$)")
    ax1.set_ylabel("Avg. nr. of correct deliveries")
    # ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    # ax1.set_yticks([0.000001, 0.0001, 0.01, 1.0, 10.0])
    # ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # # line plot on top of bar chart
    # ax2 = ax1.twinx()

    x_pos = np.arange(1, (2 * len(req_sizes)), step = 2)
    l = 0
    for cache_distance in ["1", "3"]:
        for radius in [ "S", "SF" ]:

            ax1.plot(x_pos, np.array(delivery_probs[cache_distance][radius][0]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], markersize = path_markrs_size[l], marker = path_markrs[l], label = path_labels[l])
            ax1.plot(x_pos, np.array(delivery_probs[cache_distance][radius][0]), linewidth = 1.5, color = path_markrs_colors[l], linestyle = '', markersize = path_markrs_size[l] - 3, marker = path_markrs[l], label = '')

            l = l + 1

    ax1.set_xlim(0, 10)
    ax1.set_ylim(-0.05, 2.0)

    xticks = x_pos
    xtick_labels = ['3\n', '4\n', '5\n', '10\n', '15\n']
    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels)
    # ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0])
    ax1.set_yticks([0.0, 0.5, 1.0])
    ax1.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("cache-pro-active-prefix-exclusion.pdf", bbox_inches='tight', format = 'pdf')

def plot_prefix_exclusion_table_size(data):

    event_avgs = defaultdict(defaultdict)
    path_values = defaultdict(defaultdict)
    delivery_probs = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['$|R|_{max}$ = 3', '4', '5', '10']
    event_colors = ['#000000', '#708090', '#bebebe', 'gray']
    # labels and styles for lines (path probs)
    path_labels = ['Cache pos.: 3.1, Annc. Radius: 2', 'CP: 3.1, AR: 3', 'CP: 3.3, AR: 3', 'CP: 3.3, AR: 4']
    path_colors = ['black', 'black', 'black', 'black']
    path_styles = ['-', '-', '--', '--']
    path_markrs = ['o', 'v', '^', 'o']
    path_markrs_colors = ['white', 'white', 'white', 'white']
    path_markrs_size = [10, 10, 10, 8]

    # x-axis labels
    table_sizes = []

    for key in data:

        for sub_key in data[key]:

            req_size = int(sub_key.split("-", 1)[0].split(":")[2].upper().lstrip("0"))

            if key == "path":

                table_size = int(sub_key.split("-", 1)[0].split(":")[1].upper().lstrip("0"))

                if table_size not in table_sizes:
                    table_sizes.append(table_size)


                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:4}.iteritems():

                    if k not in path_values[req_size]:
                        path_values[req_size][k] = []

                    if i not in path_probs:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[req_size][k].append(v)

                    if k == 0: 
                        delivery_probs[req_size].append(v)

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    
    # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 1

    # assumes the inter bar group space is half a bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_pos = defaultdict(list)

    show_legend = True
    for l, req_size in { 0:3, 1:4, 2:5, 3:10 }.iteritems():

        if show_legend:
            leg = event_labels[l]
        else:
            leg = None

        ax1.bar(np.arange(1, (2 * len(table_sizes)), step = 2) + (m * bar_width), np.array(delivery_probs[req_size]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
        m += 1.0

    ax1.set_ylim(0.0, 1.5)
    ax1.set_yticks([0.0, 0.5, 1.0])

    ax1.set_xlabel("Nr. of forwarding entries")
    ax1.set_ylabel("Avg. nr. of correct deliveries")
    ax1.set_xticks([1, 3, 5, 7])
    ax1.set_xticklabels(['10e5\n', '10e6\n', '10e7\n', '10e8\n']) 
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=2, loc='upper center')

    plt.savefig("cache-pro-active-prefix-exclusion-table-size.pdf", bbox_inches='tight', format = 'pdf')

def plot_global_hierarchical_entry_sizes(data):

    event_avgs = defaultdict(defaultdict)
    path_values = defaultdict(defaultdict)
    delivery_probs = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black', 'black', 'black']
    path_styles = ['-', '--', '-.', ':']
    path_markrs = ['o', 'v', '^', 'o']
    path_markrs_size = [10, 10, 10, 8]

    # x-axis labels
    entry_lengths = []

    for key in data:

        for sub_key in data[key]:

            mode = int(sub_key.split("-", 1)[1])

            if key == "events":

                entry_length = int(sub_key.split("-", 1)[0].upper().lstrip("H").lstrip("0"))
                if entry_length not in entry_lengths:
                    entry_lengths.append(entry_length)

                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:2}.iteritems():

                    if k not in event_avgs[mode]:
                        event_avgs[mode][k] = []

                    event_avgs[mode][k].append(event_probs[i])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:4}.iteritems():

                    if k not in path_values[mode]:
                        path_values[mode][k] = []

                    if i not in path_probs:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[mode][k].append(v)

    # for i in np.arange(len(path_values[1][0])):
    #     path_values[1][0][i] = path_values[1][0][i] - path_values[1][1][i]

    # for i in np.arange(len(event_avgs[1][2])):
    #     event_avgs[1][2][i] = event_avgs[1][2][i] - path_values[1][1][i] - path_values[1][2][i]

    for i in np.arange(len(path_values[1][0])):
        for s in [0, 1]:
            for mode in [0, 1]:
                delivery_probs[s].append(path_values[mode][s][i])

    for k,v in {0: 1.0, 1: 0.0}.iteritems():
        delivery_probs[k].append(v)
        delivery_probs[k].append(v)

    print("DELIVERY PROBS")
    print(delivery_probs)

    # always add the 'ideal' case: this represents the avg. number of events 
    # for a completely correct forwarding case:
    #   -# single egress iface: 6 hops
    #   -# local iface: 1 hop (the correct content source)
    #   -# no packet drops or multiple matches
    #   -# 1.0 correct deliv. prob
    #   -# 0.0 wrong deliv. prob
    for k, v in {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 0.0}.iteritems():
        for mode in [0, 1]:
            if k < 3:
                event_avgs[mode][k].append(v)
            else:
                path_values[mode][k - 3].append(v)

    entry_lengths.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    
    # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 2

    # assumes the inter bar group space is half a bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_pos = defaultdict(list)

    show_legend = True
    for mode in np.arange(2):

        if mode > 0:
            show_legend = False

        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width), np.array(event_avgs[mode][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        x_pos[mode] = np.arange(1, (2 * len(entry_lengths)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    print(interleave(x_pos[0], x_pos[1]))

    ax1.set_yscale('log')
    ax1.set_ylim(0.000001, 100000.0)

    ax1.set_xlabel("Fwd. entry size")
    ax1.set_ylabel("Avg. nr. of events.")
    # ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.set_yticks([0.000001, 0.0001, 0.01, 1.0, 10.0])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(interleave(x_pos[0], x_pos[1]), np.array(delivery_probs[l]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], markersize = path_markrs_size[l], marker = path_markrs[l], label = path_labels[l])

    ax2.set_ylim(-0.05, 2.00)
    ax2.set_ylabel("Avg. nr. of deliveries")

    xticks = [0.6, 1.0, 1.4, 2.6, 3.0, 3.4, 4.6, 5.0, 5.4, 6.6, 7.0, 7.4, 8.6, 9.0, 9.4]
    xtick_labels = ['FL', '\n1', 'R', '', '\n5', '', '', '\n10', '', '', '\n15', '', '', '\nIdeal', '']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("global-hierarchical-entry-size.pdf", bbox_inches='tight', format = 'pdf')    

def plot_global_hierarchical_bf_sizes(data):

    event_avgs = defaultdict(defaultdict)
    path_values = defaultdict(defaultdict)
    delivery_probs = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No match', 'Mult. links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe']
    # labels and styles for lines (path probs)
    path_labels = ['Corr. del.', 'Wrong del.']
    path_colors = ['black', 'black', 'black', 'black']
    path_styles = ['-', '--', '-.', ':']
    path_markrs = ['o', 'v', '^', 'o']
    path_markrs_size = [10, 10, 10, 8]

    # x-axis labels
    bf_sizes = []

    for key in data:

        for sub_key in data[key]:

            mode = int(sub_key.split("-", 1)[1][0])

            if key == "events":

                bf_size = int(sub_key.split("-", 1)[0].upper().lstrip("BF").lstrip("0"))
                if bf_size not in bf_sizes:
                    bf_sizes.append(bf_size)

                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:2}.iteritems():

                    if k not in event_avgs[mode]:
                        event_avgs[mode][k] = []

                    event_avgs[mode][k].append(event_probs[i])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                for k, i in {0:0, 1:1, 2:4}.iteritems():

                    if k not in path_values[mode]:
                        path_values[mode][k] = []

                    if i not in path_probs:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[mode][k].append(v)

    # for i in np.arange(len(path_values[1][0])):
    #     path_values[1][0][i] = path_values[1][0][i] - path_values[1][1][i] - path_values[1][2][i]

    # print("DA THINGA")
    # print(path_values[1][2])
    # print(path_values[1][1])
    # print(path_values[1][0])

    # for i in np.arange(len(event_avgs[1][2])):
    #     event_avgs[1][2][i] = event_avgs[1][2][i] - path_values[1][1][i] - path_values[1][2][i]

    for i in np.arange(len(path_values[1][0])):
        for s in [0, 1]:
            for mode in [0, 1]:
                delivery_probs[s].append(path_values[mode][s][i])

    for k,v in {0: 1.0, 1: 0.0}.iteritems():
        delivery_probs[k].append(v)
        delivery_probs[k].append(v)

    print("DELIVERY PROBS")
    print(delivery_probs)

    # always add the 'ideal' case: this represents the avg. number of events 
    # for a completely correct forwarding case:
    #   -# single egress iface: 6 hops
    #   -# local iface: 1 hop (the correct content source)
    #   -# no packet drops or multiple matches
    #   -# 1.0 correct deliv. prob
    #   -# 0.0 wrong deliv. prob
    for k, v in {0: 0.0, 1: 0.0, 2: 1.0, 3: 1.0, 4: 0.0}.iteritems():
        for mode in [0, 1]:
            if k < 3:
                event_avgs[mode][k].append(v)
            else:
                path_values[mode][k - 3].append(v)

    bf_sizes.append('Ideal')

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 4))
    ax1 = fig.add_subplot(111)
    ax1.grid(True)
    
    # we have 3 modes and 2 correctness values per mode. hence 3 groups of 
    # bars, each group with 2 bars.
    bar_group_size = len(event_labels)
    bar_group_num = 2

    # assumes the inter bar group space is half a bar. also, for n bar groups 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_pos = defaultdict(list)

    show_legend = True
    for mode in np.arange(2):

        if mode > 0:
            show_legend = False

        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(bf_sizes)), step = 2) + (m * bar_width), np.array(event_avgs[mode][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        x_pos[mode] = np.arange(1, (2 * len(bf_sizes)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    print(interleave(x_pos[0], x_pos[1]))

    ax1.set_yscale('log')
    ax1.set_ylim(0.000001, 100000.0)

    ax1.set_xlabel("BF sizes (bit)")
    ax1.set_ylabel("Avg. nr. of events.")
    # ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.set_yticks([0.000001, 0.0001, 0.01, 1.0, 10.0])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(interleave(x_pos[0], x_pos[1]), np.array(delivery_probs[l]), linewidth = 1.5, color = path_colors[l], linestyle = path_styles[l], markersize = path_markrs_size[l], marker = path_markrs[l], label = path_labels[l])

    ax2.set_ylim(-0.05, 2.0)
    ax2.set_ylabel("Avg. nr. of deliveries")

    xticks = [0.6, 1.0, 1.4, 2.6, 3.0, 3.4, 4.6, 5.0, 5.4, 6.6, 7.0, 7.4]
    xtick_labels = ['10k', '\n256', '1M', '', '\n512', '', '', '\n1024', '', '', '\nIdeal', '']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("global-hierarchical-bf-size.pdf", bbox_inches='tight', format = 'pdf')

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

def plot_pro_active_announcements(data):

    latency_avgs    = defaultdict(list)
    event_avgs      = defaultdict(defaultdict)
    path_values     = defaultdict(defaultdict)

    # labels and colors for bars (events)
    event_labels = ['Correct deliv.', 'Wrong deliv.', 'Mult. links']
    event_colors = ['lightgrey', '#708090', 'black']
    event_hatches = ['x', 'o', '*']
    # labels and styles for lines (latency probs)
    latency_labels = ['Flood', 'Random', 'Fallback', 'Ideal']
    latency_colors = ['black', 'black', 'black', 'black']
    latency_styles = ['-', '--', '-.', ':']
    latency_markrs = ['x', 'v', '^', 'o']
    latency_markrs_size = [14, 10, 10, 8]

    cache_anncs = []

    for key in data:

        for sub_key in data[key]:

            # filenames follow the format <type>.BF<bf-size>-<mode>.tsv
            # extract the mode
            mode = int(sub_key.split("-", 1)[1])

            if key == "events":

                # strip any 'BF0' prefix
                cache_annc = int(sub_key.split("-", 1)[0])

                if cache_annc not in cache_anncs:
                    cache_anncs.append(cache_annc)

                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(sub_key)
                print(event_probs)

                for i,k in {0:2, 1:2, 2:1}.iteritems():

                    if i not in event_avgs[mode]:
                        event_avgs[mode][i] = []

                    event_avgs[mode][i].append(event_probs[k])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                print(sub_key)
                print(path_probs)

                for i in [0, 1, 2]:

                    if i not in path_values[mode]:
                        path_values[mode][i] = []

                    if mode == 2 and i == 0:
                        v = path_probs[i] + path_probs[1]
                    elif i == 2:
                        v = 0.0
                    else:
                        v = path_probs[i]

                    path_values[mode][i].append(v)

                avg_latency = 0.0

                if (mode == 2):
                    latency_probs = data[key][sub_key][data[key][sub_key].STATUS <= 2]
                    latency_probs = latency_probs[latency_probs.STATUS != 1]
                else:
                    latency_probs = data[key][sub_key][data[key][sub_key].STATUS == 0]
                latency_probs = latency_probs.groupby(["LATENCY"])["PROB"].sum()
                print(latency_probs)

                # avg. latency is generally calculated as sum(latency * pron), 
                # but if the MMH_MODE is 'Flood', we simply take the min() 
                # of the keys for which prob > 0.0
                if mode == 0:

                    avg_latency = float(min(latency_probs.loc[latency_probs > 0.0].index.tolist()))

                else:

                    for lat, prob in latency_probs.iteritems():
                        avg_latency += (lat * prob)

                latency_avgs[mode].append(avg_latency)

    # print("probabilities")
    # for k in path_values:

    #     path_values[k][0] = [(a - b) for a, b in zip(path_values[k][0], path_values[k][1])]
    #     print(path_values[k][0])
    #     print(path_values[k][1])

    print(event_avgs[2][0])
    print(event_avgs[2][1])

    for mode in path_values:
        if mode == 1:
            path_values[mode][0] = [(a - b) for a,b in zip(path_values[mode][0], path_values[mode][1])]

    for mode in path_values:
        for i in np.arange(len(event_avgs[mode][2])):

            if event_avgs[mode][2][i] < 0.1:
                v = 0.1005
            else:
                v = event_avgs[mode][2][i]

            path_values[mode][2][i] = v

    # # scale the LLM avg. values by correctness
    # for mode in event_avgs:
    #     for status in event_avgs[mode]:
    #         for i in np.arange(len(event_avgs[mode][status])):

    #             if mode == 0:
    #                 event_avgs[mode][status][i] = path_values[mode][status][i]
    #             else:
    #                 event_avgs[mode][status][i] = event_avgs[mode][status][i] * path_values[mode][status][i]

    # ideal latencies are the same as with the 'Flooding' case
    for l in latency_avgs[0]:
        latency_avgs[3].append(3.0)

    print("latency")
    print(latency_avgs)

    # print("events")
    # for k in event_avgs:
    #     print(k)
    #     print(event_avgs[k])

    print("cache_anncs")
    print(cache_anncs)

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
    bar_width = 0.15

    show_legend = True
    for mode in [0, 1, 2]:

        if mode > 0:
            show_legend = False

        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(cache_anncs)), step = 2) + (m * bar_width), np.array(path_values[mode][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.1, 10.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("Cache announcement radius")
    ax1.set_ylabel("Avg. nr. of deliv.")
    # ax1.set_xticks(np.arange(1, (2 * len(cache_anncs)), step = 2), bf_sizes)
    ax1.set_yticks([0.1, 1.0, 10.0])
#    ax1.set_yticklabels(['$10^{-4}$', '10^-3', '10^-2', '10^-1', '1', '10', '', ''])

    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(latency_labels)):
        ax2.plot(np.arange(1, (2 * len(cache_anncs)), step = 2), np.array(latency_avgs[l]), linewidth = 1.5, color = latency_colors[l], linestyle = latency_styles[l], marker = latency_markrs[l], markersize = latency_markrs_size[l], label = latency_labels[l])

    ax2.set_ylim(0.0, 18.0)
    ax2.set_ylabel("Avg. latency")
    xticks = [0.4, 1.0, 1.6, 2.4, 3.0, 3.6, 4.4, 5.0, 5.6]
    # xtick_labels = ['FL', 'R\n256', 'FB', 'FL', 'R\n512', 'FB', 'FL', 'R\n1024', 'FB']
    xtick_labels = ['FL', 'R\n1', 'FB', '', '\n2', '', '', '\n3', '']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0])
    # ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("cache-pro-active-annc.pdf", bbox_inches='tight', format = 'pdf')

def plot_opportunistic_on_path(data):

    latency_avgs    = defaultdict(list)
    event_avgs      = defaultdict(defaultdict)
    path_values     = defaultdict(defaultdict)

    # labels and colors for bars (events)
    event_labels = ['Correct deliv.', 'Wrong deliv.', 'No match']
    event_colors = ['lightgrey', '#708090', 'black']
    event_hatches = ['x', 'o', '*']
    # labels and styles for lines (latency probs)
    latency_labels = ['Random', 'Fallback', 'Ideal']
    latency_colors = ['black', 'black', 'black']
    latency_styles = ['-', '--', '-.']
    latency_markrs = ['x', 'v', 'o']
    latency_markrs_size = [14, 10, 8]

    cache_distances = []
    probs_sum = defaultdict(list)

    for key in data:
        for sub_key in data[key]:

            # filenames follow the format: 
            # <type>.<cache-annc>:<mmh-mode>.tsv

            # extract mult. match handling mode
            mmh_mode = sub_key.split(":", 1)[1].upper()

            if key == "events":

                # extract the cache dist and bf size
                cache_distance = sub_key.split(":", 1)[0].upper()

                if cache_distance not in cache_distances:
                    cache_distances.append(cache_distance)

                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(sub_key)
                print(event_probs)

                for i, e in {0:2, 1:2, 2:0}.iteritems():

                    if i not in event_avgs[mmh_mode]:
                        event_avgs[mmh_mode][i] = []
                    
                    event_avgs[mmh_mode][i].append(event_probs[e])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                print(sub_key)
                print(path_probs)

                p = 0.0
                for i in [0, 1]:

                    v = path_probs[i]

                    if i not in path_values[mmh_mode]:
                        path_values[mmh_mode][i] = []

                    path_values[mmh_mode][i].append(v)
                    p += v

                probs_sum[mmh_mode].append(p)

                avg_latency = 0.0
                latency_probs = data[key][sub_key].groupby(by = ["LATENCY"])["PROB"].sum()

                for lat, prob in latency_probs.iteritems():
                    avg_latency += (lat * prob)

                latency_avgs[mmh_mode].append(avg_latency)

    print("probs")
    for k in path_values:
        for i in [0, 1]:
            print(k)
            print(path_values[k][i])
            print(probs_sum[k])
            path_values[k][i] = [ a / b for a, b in zip(path_values[k][i], probs_sum[k]) ]
            print(path_values[k][i])

    # scale the LLM avg. values by correctness
    for mmh_mode in event_avgs:
        for status in event_avgs[mmh_mode]:
            for i in np.arange(len(event_avgs[mmh_mode][status])):

                if status < 2:
                    v = event_avgs[mmh_mode][status][i] * path_values[mmh_mode][status][i]
                else:
                    v = event_avgs[mmh_mode][status][i]

                # ok, this is cheating a bit: the idea is to show a bit 
                # of a bar (representing 0.0), just for the reader to know that this 
                # isn't an error
                if v < 0.00001:
                    v = 0.0000105

                event_avgs[mmh_mode][status][i] = v

    # ideal latencies (added manually)
    latency_avgs["I"].append(2.0)
    latency_avgs["I"].append(3.0)
    latency_avgs["I"].append(3.0)

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
    for mmh_mode in event_avgs:

        # print the bars for each |L|, for both cache distances
        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(cache_distances)), step = 2) + (m * bar_width), np.array(event_avgs[mmh_mode][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        # show legend only once
        show_legend = False

        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.00001, 1000.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("Cache distance")
    ax1.set_ylabel("Avg. nr. of deliv.")
    # ax1.set_xticks(np.arange(1, (2 * len(cache_distances)), step = 2), cache_distances)
    # ax1.set_yticks(np.arange(0.0, 1.75, step = 0.25))
    # ax1.set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0', '1.25', ''])
    ax1.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for k, l in {0:"R", 1:"F", 2:"I"}.iteritems():

        ax2.plot(np.arange(1, (2 * len(cache_distances)), step = 2), np.array(latency_avgs[l]), linewidth = 1.5, color = latency_colors[k], linestyle = latency_styles[k], marker = latency_markrs[k], markersize = latency_markrs_size[k], label = latency_labels[k])

    ax2.set_ylim(0.0, 14.0)
    ax2.set_ylabel("Avg. latency")

    xticks = [0.6, 1.0, 1.4, 2.6, 3.0, 3.6, 4.6, 5.0, 5.4]
    xtick_labels = ['R', '\n2', 'FB', 'R', '\n3', 'FB', 'R', '\n3 (512 BF)', 'FB']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    # ax2.set_xlim(x_min, x_max)

    ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
    # ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("cache-opportunistic-on-path.pdf", bbox_inches='tight', format = 'pdf')

def plot_opportunistic_off_path(data):

    latency_avgs    = defaultdict(list)
    event_avgs      = defaultdict(defaultdict)
    path_values     = defaultdict(defaultdict)

    # labels and colors for bars (events)
    event_labels = ['Correct deliv.', 'Wrong deliv.']
    event_colors = ['lightgrey', '#708090']
    event_hatches = ['x', 'o', '*']
    # labels and styles for lines (latency probs)
    latency_labels = ['Random', 'Fallback', 'Ideal']
    latency_colors = ['black', 'black', 'black']
    latency_styles = ['-', '--', '-.']
    latency_markrs = ['x', 'v', 'o']
    latency_markrs_size = [14, 10, 8]

    cache_anncs = []
    probs_sum = defaultdict(list)

    for key in data:
        for sub_key in data[key]:

            # filenames follow the format: 
            # <type>.<cache-annc>:<mmh-mode>.tsv

            # extract mult. match handling mode
            mmh_mode = sub_key.split(":", 1)[1].upper()

            if key == "events":

                # extract the cache dist and bf size
                cache_annc = int(sub_key.split(":", 1)[0].upper())

                if cache_annc not in cache_anncs:
                    cache_anncs.append(cache_annc)

                # group by event number and sum the probabilities to get the 
                # expected value (i.e. avg. nr. of events by type)
                event_probs = data[key][sub_key].groupby(by = ["EVENT"])["PROB"].sum()

                print(sub_key)
                print(event_probs)

                for i, e in {0:2, 1:2}.iteritems():

                    if i not in event_avgs[mmh_mode]:
                        event_avgs[mmh_mode][i] = []
                    
                    event_avgs[mmh_mode][i].append(event_probs[e])

            elif key == "path":

                path_probs = data[key][sub_key].groupby(by = ["STATUS"])["PROB"].sum()

                print(sub_key)
                print(path_probs)

                p = 0.0
                for i in [0, 1]:

                    v = path_probs[i]

                    if i not in path_values[mmh_mode]:
                        path_values[mmh_mode][i] = []

                    path_values[mmh_mode][i].append(v)
                    p += v

                probs_sum[mmh_mode].append(p)

                avg_latency = 0.0
                latency_probs = data[key][sub_key].groupby(by = ["LATENCY"])["PROB"].sum()

                for lat, prob in latency_probs.iteritems():
                    avg_latency += (lat * prob)

                latency_avgs[mmh_mode].append(avg_latency)

    print("probs")
    for k in path_values:
        for i in [0, 1]:
            print(k)
            print(path_values[k][i])
            print(probs_sum[k])
            path_values[k][i] = [ a / b for a, b in zip(path_values[k][i], probs_sum[k]) ]
            print(path_values[k][i])

    # scale the LLM avg. values by correctness
    for mmh_mode in event_avgs:
        for status in event_avgs[mmh_mode]:
            for i in np.arange(len(event_avgs[mmh_mode][status])):

                v = event_avgs[mmh_mode][status][i] * path_values[mmh_mode][status][i]

                # # ok, this is cheating a bit: the idea is to show a bit 
                # # of a bar (representing 0.0), just for the reader to know that this 
                # # isn't an error
                # if v < 0.005:
                #     v = 0.005

                event_avgs[mmh_mode][status][i] = v

    # ideal latencies (added manually)
    latency_avgs["I"].append(3.0)
    latency_avgs["I"].append(3.0)
    latency_avgs["I"].append(3.0)

    print("latency")
    print(latency_avgs)

    print("events")
    for k in event_avgs:
        print(k)
        print(event_avgs[k])

    print("cache anncs")
    print(cache_anncs)

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

    bar_group_offset = np.arange(1, (2 * len(cache_anncs)), step = 2)[0] + (m * bar_width)
    bar_group_width = (bar_group_num * bar_group_size + 1) * bar_width

    x_max = np.arange(1, (2 * len(cache_anncs)), step = 2)[0] + (m * bar_width)
    x_max += np.arange(1, (2 * len(cache_anncs)), step = 2)[-1] + ((bar_group_num * bar_group_size + 1) / 2.0) * bar_width

    print(bar_group_size)
    print(bar_group_num)
    print(bar_group_offset)
    print(bar_group_width)
    print("[%d, %d]" % (x_min, x_max))

    show_legend = True
    for mmh_mode in event_avgs:

        # print the bars for each |L|, for both cache distances
        for l in np.arange(len(event_labels)):

            if show_legend:
                leg = event_labels[l]
            else:
                leg = None

            ax1.bar(np.arange(1, (2 * len(cache_anncs)), step = 2) + (m * bar_width), np.array(event_avgs[mmh_mode][l]), color = event_colors[l], linewidth = 1.5, alpha = 0.75, width = bar_width, label = leg)
            m += 1.0

        # show legend only once
        show_legend = False

        m += 1.0

    ax1.set_yscale('log')
    ax1.set_ylim(0.00001, 1000.0)

    # ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("Cache announcement radius")
    ax1.set_ylabel("Avg. nr. of deliv.")
    # ax1.set_xticks(np.arange(1, (2 * len(cache_distances)), step = 2), cache_distances)
    # ax1.set_yticks(np.arange(0.0, 1.75, step = 0.25))
    # ax1.set_yticklabels(['0.0', '0.25', '0.5', '0.75', '1.0', '1.25', ''])
    ax1.set_yticks([0.00001, 0.0001, 0.001, 0.01, 0.1, 1.0])
    ax1.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for k, l in {0:"R", 1:"F", 2:"I"}.iteritems():

        ax2.plot(np.arange(1, (2 * len(cache_anncs)), step = 2), np.array(latency_avgs[l]), linewidth = 1.5, color = latency_colors[k], linestyle = latency_styles[k], marker = latency_markrs[k], markersize = latency_markrs_size[k], label = latency_labels[k])

    ax2.set_ylim(0.0, 14.0)
    ax2.set_ylabel("Avg. latency")

    xticks = [0.7, 1.0, 1.3, 2.7, 3.0, 3.3, 4.7, 5.0, 5.3]
    xtick_labels = ['R', '\n0', 'FB', 'R', '\n1', 'FB', 'R', '\n2', 'FB']
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)
    # ax2.set_xlim(x_min, x_max)

    ax2.set_yticks([0.0, 2.0, 4.0, 6.0, 8.0, 10.0, 12.0, 14.0])
    # ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=FONTSIZE_LEGEND, ncol=1, loc='upper right')

    plt.savefig("cache-opportunistic-off-path.pdf", bbox_inches='tight', format = 'pdf')

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

    if args.case == 'global':

        if args.subcase in ['flood', 'random']:
            plot_global(data, args.subcase)

        elif args.subcase.split("-", 1)[0] == 'topologies':

            if args.subcase.split("-", 1)[1] == 'entry':
                plot_global_topologies_entry_sizes(data)
            elif args.subcase.split("-", 1)[1] == 'bf-sizes':
                plot_global_topologies_bf_sizes(data)
            else:
                sys.stderr.write("""%s: [ERROR] please supply a valid subcase ('topologies-entry' or 'topologies-bf-sizes').\n""" % sys.argv[0]) 
                parser.print_help()
                sys.exit(1)

        elif args.subcase.split("-", 1)[0] == 'hierarchical':

            if args.subcase.split("-", 1)[1] == 'entry':
                plot_global_hierarchical_entry_sizes(data)
            else:
                plot_global_hierarchical_bf_sizes(data)

        elif args.subcase == 'peering':
            plot_global_peering(data)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase ('flood', 'random', 'hierarchical' or 'peering').\n""" % sys.argv[0]) 
            parser.print_help()
            sys.exit(1)            

    elif args.case == 'bf-sizes':
        plot_bf_sizes(data)
    elif args.case == 'req-sizes':
        plot_req_sizes(data)

    elif args.case == 'cache-pro-active':

        if args.subcase == 'anncs':
            plot_pro_active_announcements(data)
        elif args.subcase == 'bf-sizes':
            plot_pro_active_latencies(data)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase ('anncs' or 'bf-sizes').\n""" % sys.argv[0]) 
            parser.print_help()
            sys.exit(1)

    elif args.case == 'cache-opportunistic':

        if args.subcase == 'on-path':
            plot_opportunistic_on_path(data)
        elif args.subcase == 'off-path':
            plot_opportunistic_off_path(data)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase ('on-path' or 'off-path').\n""" % sys.argv[0]) 
            parser.print_help()
            sys.exit(1)

    elif args.case == 'prefix-exclusion':

        if args.subcase == 'req-sizes':
            plot_prefix_exclusion(data)
        elif args.subcase == 'table-sizes':
            plot_prefix_exclusion_table_size(data)
        else:
            sys.stderr.write("""%s: [ERROR] please supply a valid subcase ('req_sizes').\n""" % sys.argv[0]) 
            parser.print_help()
            sys.exit(1)

    else:
        sys.stderr.write("""%s: [ERROR] please supply a valid case\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)