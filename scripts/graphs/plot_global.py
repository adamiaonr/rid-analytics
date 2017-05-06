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

def get_path_info(test_file, file_label):

    # extract info from .test file
    test_run = et.parse(test_file)
    test_run_root = test_run.getroot()

    # use test_id to get additional information about the 
    # test from the .test file
    test_id = ""
    if 'R' in file_label:
        test_id = '-'.join(file_label.split("-")[:8])
    else:
        test_id = '-'.join(file_label.split("-")[:7])
    # print("get_path_length() : test_id = %s" % (test_id))

    for test in test_run.findall('test'):
        # test id to use as prefix to output labels
        if test.get('id') != test_id:
            continue

        for path in test.find('paths').findall('path'):
            return (len(path.text.split(',')) - 1), float(path.get('avg_outdegree'))

def plot_bf_sizes_fallbacks(data, test_dir):

    # forwarding efficiency values
    fwd_events = defaultdict()
    fwd_efficiency = defaultdict()
    # latency averages
    avg_delivery_latencies = defaultdict()
    # delivery probs
    avg_delivery_probs = defaultdict()

    # labels and colors for bars (topologies)
    topology_keys = ['1221', '4755', '7018']
    topology_colors = ['#000000', '#708090', '#bebebe']
    radius_keys = ['1', '2']

    # modes : '0' flood, '2' fallback
    modes = []

    # path info
    path_lengths = defaultdict()
    avg_outdegrees = defaultdict()

    # populate the fwd_* dicts w/ data obtained from .tsv files
    got_path_info = []

    for file_type in data:
        for file_label in data[file_type]:

            # collect topology and bf-size keys
            topology_key = file_label.split("-")[0]
            radius_key = file_label.split("-")[7].lstrip("R")

            if topology_key not in got_path_info:
                path_lengths[topology_key], avg_outdegrees[topology_key] = get_path_info(os.path.join(test_dir, ("%s.test" % (topology_key))), file_label)
                got_path_info.append(topology_key)

            # collect the mode as int
            mode = int(file_label.split("-")[5])
            if mode not in modes:
                modes.append(mode)

            if file_type == "events":

                # gather the total probs by event
                event_probs = data[file_type][file_label].groupby(by = ["EVENT"])["PROB"].sum()

                # gather the probabilities for MLM and SLM events for 
                # forwarding efficiency (these events are mutually exclusive)
                if topology_key not in fwd_efficiency:
                    fwd_efficiency[topology_key] = defaultdict()
                    fwd_events[topology_key] = defaultdict()

                if radius_key not in fwd_efficiency[topology_key]:
                    fwd_efficiency[topology_key][radius_key] = []
                    fwd_events[topology_key][radius_key] = []

                events = event_probs[EVENT_MLM] + event_probs[EVENT_SLM]

                fwd_events[topology_key][radius_key].append(events)
                fwd_efficiency[topology_key][radius_key].append(float(path_lengths[topology_key]) / float(events))

            if file_type == "path":

                # this will save the average delivery latency value
                avg_delivery_latency = 0.0
                # keep track of the delivery probability
                delivery_prob = 0.0
                delivery_probs = data[file_type][file_label].groupby(by = ["STATUS"])["PROB"].sum()

                if topology_key not in avg_delivery_latencies:
                    avg_delivery_latencies[topology_key] = defaultdict()
                    avg_delivery_probs[topology_key] = defaultdict()

                if radius_key not in avg_delivery_latencies[topology_key]:
                    avg_delivery_latencies[topology_key][radius_key] = []
                    avg_delivery_probs[topology_key][radius_key] = []

                avg_delivery_probs[topology_key][radius_key].append((delivery_probs[0], delivery_probs[1]))

                # order dataframe values by STATUS and LATENCY. the objective is 
                # to only look at the STATUS = 0 (correct delivery) and keep 
                # multiplying LATENCY * PROB until we reach a delivery prob of 1.0
                delivery_latency = data[file_type][file_label].sort(["STATUS","LATENCY"])
                
                for index, row in delivery_latency.iterrows():

                    print("%s.%s.%d : %d, %d, %.2f" % (topology_key, radius_key, mode, row['STATUS'], row['LATENCY'], row['PROB']))
                    if (delivery_prob + row['PROB']) >= 1.0:
                        avg_delivery_latency += ((1.0 - delivery_prob) * row['LATENCY'])
                        break

                    delivery_prob += row['PROB']
                    avg_delivery_latency += (row['PROB'] * row['LATENCY'])

                avg_delivery_latencies[topology_key][radius_key].append(avg_delivery_latency)

    print("fwd_efficiency : %s" % (fwd_efficiency))
    print("fwd_events :")
    for t in fwd_events:
        for b in fwd_events[t]:
            print("%s.%s = %s" % (t, b, str(fwd_events[t][b])))

    print("avg_delivery_latencies :")
    for t in avg_delivery_latencies:
        for b in avg_delivery_latencies[t]:
            print("%s.%s = %s" % (t, b, str(avg_delivery_latencies[t][b])))

    print("avg_delivery_probs :")
    for t in avg_delivery_probs:
        for b in avg_delivery_probs[t]:
            print("%s.%s = %s" % (t, b, str(avg_delivery_probs[t][b])))

    for t in path_lengths:
        print("%s : %d" % (t, path_lengths[t]))

    # create labels for the bars
    topology_labels = {'1221': 'Telstra (1221)', '3257': 'Tiscali (3257)', '4755': 'VSNL (4755', '7018': 'AT&T (7018)'}

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True)

    bar_group_size = len(topology_keys)
    bar_group_num = 2

    # assumes the inter bar group space is half a bar. also, for n groups of bars 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.15

    x_pos = defaultdict(list)
    xx_pos = defaultdict(list)

    show_legend = True
    i = 0
    for radius_key in radius_keys:

        if radius_key != '1':
            show_legend = False

        for t, topology_key in enumerate(topology_keys):

            if show_legend:
                leg = topology_labels[topology_key]
            else:
                leg = None

            xx_pos[i] = (np.arange(1, (2 * len(modes)), step = 2) + (m * bar_width) + (bar_width / 2.0))
            i += 1
            ax1.bar(np.arange(1, (2 * len(modes)), step = 2) + (m * bar_width), np.array(avg_delivery_latencies[topology_key][radius_key]), color = topology_colors[t], linewidth = 1.5, alpha = 0.55, width = bar_width, label = leg)
            m += 1.0

        x_pos[radius_key] = np.arange(1, (2 * len(modes)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    ax1.set_xlabel("Annc. radius\nMultiple match res. mode")
    ax1.set_ylabel("Avg. latency (# of hops)")
    # ax1.set_yscale('log')
    ax1.set_ylim(0, 12)
    ax1.set_yticks(np.arange(0, 8 + 1, step = 2))

    xxticks = x_pos['1'] + ((x_pos['2'][0] - x_pos['1'][0]) / 2.0)
    xticks = interleave_n(x_pos['1'], xxticks, x_pos['2'])
    # a convoluted way to set the x-axis labels?
    xtick_labels = []
    for mode_str in ['Flood', 'Fallback']:
        xtick_labels.append("%s" % ('1'))
        xtick_labels.append("\n%s" % (mode_str))
        xtick_labels.append("%s" % ('2'))

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # fwd_efficiency as a line plot on top of bar chart
    ax2 = ax1.twinx()
    ax2.yaxis.grid(False)
    # a convoluted way to transform a dictionary into a list, 
    # ready to use in plot(). the outermost guide of the cycle 
    # is the entry size |f|
    avg_correct_delivery_probs_array = []
    avg_incorrect_delivery_probs_array = []

    for f, el in enumerate(modes):
        for radius_key in radius_keys:
            for topology_key in topology_keys:
                avg_correct_delivery_probs_array.append(int(avg_delivery_probs[topology_key][radius_key][f][0]))
                avg_incorrect_delivery_probs_array.append(int(avg_delivery_probs[topology_key][radius_key][f][1]))

    # another convoluted way to create the x-axis for the fwd efficiency graph
    # FIXME: interleave_n() accepts a list of lists
    print(xx_pos)
    xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5])

    ax2.plot(xx_pos, avg_correct_delivery_probs_array, linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = '^', label = 'Corr. del.')
    ax2.plot(xx_pos, avg_incorrect_delivery_probs_array, linewidth = 1.5, color = 'black', linestyle = '--', markersize = 5, marker = 'v', label = 'Incorr. del.')
    # ax2.axhspan(0, 100, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    ax2.set_xlim(xx_pos[0] - (3 * bar_width), xx_pos[-1] + (3 * bar_width))
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)

    ax2.set_ylim(-20, 40)
    ax2.set_yticks([0, 10, 20])
    ax2.set_ylabel("Avg. # of deliveries")

    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("cdn-fallbacks.pdf", bbox_inches='tight', format = 'pdf')

def plot_request_sizes(data, test_dir):

    fwd_events = defaultdict()
    fwd_efficiency = defaultdict()

    # labels and colors for bars (topologies)
    topology_keys = ['1221', '4755', '7018']
    topology_colors = ['#000000', '#708090', '#bebebe']
    bf_keys = ['192', '256']

    # x-axis labels (|f| values)
    req_sizes = []

    # path info
    path_lengths = defaultdict()
    avg_outdegrees = defaultdict()

    # populate the fwd_* dicts w/ data obtained from .tsv files
    got_path_info = []

    for file_type in data:
        for file_label in data[file_type]:

            # collect topology and bf-size keys
            topology_key = file_label.split("-")[0]
            bf_key = file_label.split("-")[1]

            if topology_key not in got_path_info:
                path_lengths[topology_key], avg_outdegrees[topology_key] = get_path_info(os.path.join(test_dir, ("%s.test" % (topology_key))), file_label)
                got_path_info.append(topology_key)

            req_size = int(file_label.split("-")[2])
            if req_size not in req_sizes:
                # print("scanning |f| = %d" % (req_size))
                req_sizes.append(req_size)

            if file_type == "events":

                # gather the total probs by event
                event_probs = data[file_type][file_label].groupby(by = ["EVENT"])["PROB"].sum()

                # gather the probabilities for MLM and SLM events for 
                # forwarding efficiency (these events are mutually exclusive)
                if topology_key not in fwd_efficiency:
                    fwd_efficiency[topology_key] = defaultdict()
                    fwd_events[topology_key] = defaultdict()

                if bf_key not in fwd_efficiency[topology_key]:
                    fwd_efficiency[topology_key][bf_key] = []
                    fwd_events[topology_key][bf_key] = []

                events = event_probs[EVENT_MLM] + event_probs[EVENT_SLM]

                fwd_events[topology_key][bf_key].append(events)
                fwd_efficiency[topology_key][bf_key].append(float(path_lengths[topology_key]) / float(events))

    print("fwd_efficiency : %s" % (fwd_efficiency))
    print("fwd_events :")
    for t in fwd_events:
        for b in fwd_events[t]:
            print("%s.%s = %s" % (t, b, str(fwd_events[t][b])))

    for t in path_lengths:
        print("%s : %d" % (t, path_lengths[t]))

    # create labels for the bars
    topology_labels = {'1221': 'Telstra (1221)', '3257': 'Tiscali (3257)', '4755': 'VSNL (4755', '7018': 'AT&T (7018)'}

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True)

    # each |f| slot will have 3 groups of bars, w/ 
    # 3 bars each. there will be a total of 
    # 18 bars in the graph
    bar_group_size = len(topology_keys)
    bar_group_num = 3

    # assumes the inter bar group space is half a bar. also, for n groups of bars 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.15

    x_pos = defaultdict(list)
    xx_pos = defaultdict(list)

    show_legend = True
    i = 0
    for bf_key in bf_keys:

        if bf_key != '192':
            show_legend = False

        for t, topology_key in enumerate(topology_keys):

            if show_legend:
                leg = topology_labels[topology_key]
            else:
                leg = None

            xx_pos[i] = (np.arange(1, (2 * len(req_sizes)), step = 2) + (m * bar_width) + (bar_width / 2.0))
            i += 1
            ax1.bar(np.arange(1, (2 * len(req_sizes)), step = 2) + (m * bar_width), np.array(fwd_events[topology_key][bf_key]), color = topology_colors[t], linewidth = 1.5, alpha = 0.55, width = bar_width, label = leg)
            m += 1.0

        x_pos[bf_key] = np.arange(1, (2 * len(req_sizes)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    ax1.set_xlabel("BF sizes\nRequest size")
    ax1.set_ylabel("Avg. # of used links")
    # ax1.set_yscale('log')
    ax1.set_ylim(0, 50)
    ax1.set_yticks([0, 10, 20, 30])

    xxticks = x_pos['192'] + ((x_pos['256'][0] - x_pos['192'][0]) / 2.0)
    xticks = interleave_n(x_pos['192'], xxticks, x_pos['256'])
    # a convoluted way to set the x-axis labels?
    xtick_labels = []
    for req_size in req_sizes:
        xtick_labels.append("%s" % ('192'))
        xtick_labels.append("\n%d" % (req_size))
        xtick_labels.append("%s" % ('256'))

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # fwd_efficiency as a line plot on top of bar chart
    ax2 = ax1.twinx()
    ax2.yaxis.grid(False)
    # a convoluted way to transform a dictionary into a list, 
    # ready to use in plot(). the outermost guide of the cycle 
    # is the entry size |f|
    fwd_efficiency_array = []
    for f, el in enumerate(req_sizes):
        for bf_key in bf_keys:
            for topology_key in topology_keys:

                fwd_efficiency_array.append(int(fwd_efficiency[topology_key][bf_key][f] * 100.0))

    # another convoluted way to create the x-axis for the fwd efficiency graph
    # FIXME: interleave_n() accepts a list of lists
    print(xx_pos)
    xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5])

    ax2.plot(xx_pos, fwd_efficiency_array, linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = 'Fwd. effic.')
    # ax2.axhspan(0, 100, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    ax2.set_xlim(xx_pos[0] - (3 * bar_width), xx_pos[-1] + (3 * bar_width))
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)

    ax2.set_ylim(-50, 200)
    ax2.set_yticks([0, 50, 100])
    ax2.set_ylabel("Fwd. efficiency (%)")

    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("global-req-size-efficiency.pdf", bbox_inches='tight', format = 'pdf')

def plot_bf_sizes(data, test_dir):

    # we plot 2 things : forwarding efficiency and avg. nr. of deliveries. the 
    # objective is to plot a bar chart for fwd efficiency values, together with 
    # a line chart for avg. nr. of deliveries. the look of the graph 
    # should be as below:
    #
    #  #: length 1, %: length 3, &: length 6
    #
    #    *                    |- 6 (avg. deliveries)
    #  *- -*                  |- 5
    #       \   /*-*\    *-*  |- 4
    #    %   --*     \-*/     |- 3
    #  # %       %            |- 2
    #  # % &   # % 	     %    |- 1
    #  # % &   # % &   # % &  |- 0
    # ----------------------- 
    #  |f|=1   |f|=5   |f|=10
    #

    # the efficiency values are saved in a dict, indexed by the length of the 
    # path (1, 3 and 6 hops). each key points to a list of 3 values, for the 3 
    # diff. entry sizes : 1, 5 and 10. 
    fwd_events = defaultdict()
    fwd_efficiency = defaultdict()

    # labels and colors for bars (topologies)
    topology_keys = ['1221', '4755']
    topology_colors = ['#000000', '#708090', '#bebebe']

    # x-axis labels (|f| values)
    entry_lengths = [1, 5]

    # path info
    path_lengths = defaultdict()
    avg_outdegrees = defaultdict()

    # nr. of samples
    num_samples = defaultdict()

    # populate the fwd_* dicts w/ data obtained from .tsv files
    got_path_info = []

    for file_type in data:
        for file_label in data[file_type]:

            # collect topology and bf-size keys
            topology_key = file_label.split("-")[0]
            bf_key = file_label.split("-")[1]

            if topology_key not in got_path_info:
                path_lengths[topology_key], avg_outdegrees[topology_key] = get_path_info(os.path.join(test_dir, ("%s.test" % (topology_key))), file_label)
                got_path_info.append(topology_key)

            entry_length = int(file_label.split("-")[3])
            # if entry_length not in entry_lengths:
            #     # print("scanning |f| = %d" % (entry_length))
            #     entry_lengths.append(entry_length)

            if file_type == "events":

                # gather the total probs by event
                event_probs = data[file_type][file_label].groupby(by = ["EVENT"])["PROB"].sum()

                # gather the probabilities for MLM and SLM events for 
                # forwarding efficiency (these events are mutually exclusive)
                if topology_key not in fwd_efficiency:
                    fwd_efficiency[topology_key] = defaultdict()
                    fwd_events[topology_key] = defaultdict()
                    num_samples[topology_key] = defaultdict()

                if bf_key not in fwd_efficiency[topology_key]:
                    fwd_efficiency[topology_key][bf_key] = [0.0, 0.0]
                    fwd_events[topology_key][bf_key] = [0.0, 0.0]
                    num_samples[topology_key][bf_key] = [0, 0]

                events = event_probs[EVENT_MLM] + event_probs[EVENT_SLM]

                # if len(fwd_events[topology_key][bf_key]) < 2:
                #     fwd_events[topology_key][bf_key].append(events)
                #     fwd_efficiency[topology_key][bf_key].append(float(path_lengths[topology_key]) / float(events))
                #     num_samples[topology_key][bf_key].append(1)
                # else:
                fwd_events[topology_key][bf_key][entry_lengths.index(entry_length)] += events
                print("fwd_efficiency %s:%s:%d : %f" % (topology_key, bf_key, entry_length, (float(path_lengths[topology_key]) / float(events))))
                fwd_efficiency[topology_key][bf_key][entry_lengths.index(entry_length)] += (float(path_lengths[topology_key]) / float(events))
                num_samples[topology_key][bf_key][entry_lengths.index(entry_length)] += 1

                print("num_samples %s:%s:%d : %s" % (topology_key, bf_key, entry_length, str(num_samples[topology_key][bf_key])))


    print("fwd_efficiency :")
    for t in fwd_efficiency:
        for b in fwd_efficiency[t]:
            print("%s.%s = %s" % (t, b, str(fwd_efficiency[t][b])))

    print("fwd_events :")
    for t in fwd_events:
        for b in fwd_events[t]:
            print("%s.%s = %s" % (t, b, str(fwd_events[t][b])))

    print("num_samples :")
    for t in num_samples:
        for b in num_samples[t]:
            print("%s.%s = %s" % (t, b, str(num_samples[t][b])))

    print("path_lengths :")
    for t in path_lengths:
        print("%s : %d" % (t, path_lengths[t]))

    # create labels for the bars
    topology_labels = {'1221': 'Telstra (1221)', '3257': 'Tiscali (3257)', '4755': 'VSNL (4755)', '7018': 'AT&T (7018)'}

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True)

    # each |f| slot will have 3 groups of bars, w/ 
    # 3 bars each. there will be a total of 
    # 18 bars in the graph
    bar_group_size = len(topology_keys)
    bar_group_num = 3

    # assumes the inter bar group space is half a bar. also, for n groups of bars 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 2.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.15

    x_pos = defaultdict(list)
    xx_pos = defaultdict(list)

    show_legend = True
    i = 0
    for bf_key in ['192', '256', '512']:

        if bf_key != '192':
            show_legend = False

        for t, topology_key in enumerate(topology_keys):

            if show_legend:
                leg = topology_labels[topology_key]
            else:
                leg = None

            xx_pos[i] = (np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width) + (bar_width / 2.0))
            i += 1
            ax1.bar(
                np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width), 
                np.array(fwd_events[topology_key][bf_key]) / np.array(num_samples[topology_key][bf_key]), 
                color = topology_colors[t], linewidth = 1.5, alpha = 0.55, width = bar_width, label = leg)
            m += 1.0

        x_pos[bf_key] = np.arange(1, (2 * len(entry_lengths)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    ax1.set_xlabel("BF sizes\nFwd. entry size")
    ax1.set_ylabel("Avg. # of used links")
    # ax1.set_yscale('log')
    ax1.set_ylim(0, 200)
    ax1.set_yticks([0, 20, 40, 60, 80, 100, 120, 140])

    xticks = interleave_n(x_pos['192'], x_pos['256'], x_pos['512'])
    # a convoluted way to set the x-axis labels?
    xtick_labels = []
    for entry_length in entry_lengths:
        for bf_key in ['192', '256', '512']:
            if bf_key == '256':
                xtick_labels.append("%d\n%d" % (int(bf_key), entry_length))
            else:
                xtick_labels.append("%d" % (int(bf_key)))

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # fwd_efficiency as a line plot on top of bar chart
    ax2 = ax1.twinx()
    ax2.yaxis.grid(False)
    # a convoluted way to transform a dictionary into a list, 
    # ready to use in plot(). the outermost guide of the cycle 
    # is the entry size |f|
    fwd_efficiency_array = []
    for f, el in enumerate(entry_lengths):
        for bf_key in ['192', '256', '512']:
            for topology_key in topology_keys:

                fwd_efficiency_array.append(int((fwd_efficiency[topology_key][bf_key][f] / num_samples[topology_key][bf_key][f]) * 100.0))

    # another convoluted way to create the x-axis for the fwd efficiency graph
    # FIXME: interleave_n() accepts a list of lists
    print(xx_pos)
    xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5])
    # xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5], xx_pos[6], xx_pos[7], xx_pos[8])

    ax2.plot(xx_pos, fwd_efficiency_array, linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = 'Fwd. effic.')
    # ax2.axhspan(0, 100, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    ax2.set_xlim(xx_pos[0] - (2 * bar_width), xx_pos[-1] + (2 * bar_width))
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)

    ax2.set_ylim(-50, 200)
    ax2.set_yticks([0, 50, 100])
    ax2.set_ylabel("Fwd. efficiency (%)")

    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("global-bf-sizes.pdf", bbox_inches='tight', format = 'pdf') 