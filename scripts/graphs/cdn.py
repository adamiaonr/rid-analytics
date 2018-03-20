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

def plot_fallbacks(data, test_dir):

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
    bf_keys = ['192', '256']

    # modes : '1' relf, '2' fallback
    modes = [1, 2]

    # annc. size
    annc_sizes = [2, 5]

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
            if bf_key == '192':
                continue

            # if topology_key not in got_path_info:
            #     path_lengths[topology_key], avg_outdegrees[topology_key] = get_path_info(os.path.join(test_dir, ("%s.test" % (topology_key))), file_label)
            #     got_path_info.append(topology_key)

            # collect the mode as int
            mode = int(file_label.split("-")[5])
            # if mode not in modes:
            #     modes.append(mode)

            # 14-2, 37-2

            annc_size = 0

            if not ('tp2' in file_label or 'tp5' in file_label):
                continue

            if file_label.split("-")[7] != 'tp5':

                annc_size = 2

                # src_dst_pair = (int(file_label.split("-")[7]), int(file_label.split("-")[8]))
                # if topology_key == '1221' and (src_dst_pair == (14,2) or src_dst_pair == (37,2)):
                #     continue

                # if topology_key == '7018' and (src_dst_pair == (86,10)):
                #     continue

            else:

                annc_size = 5


            if file_type == "events":

                # gather the total probs by event
                event_probs = data[file_type][file_label].groupby(by = ["EVENT"])["PROB"].sum()

                # gather the probabilities for MLM and SLM events for 
                # forwarding efficiency (these events are mutually exclusive)
                if topology_key not in fwd_efficiency:
                    fwd_efficiency[topology_key] = defaultdict()
                    fwd_events[topology_key] = defaultdict()
                    num_samples[topology_key] = defaultdict()

                if mode not in fwd_efficiency[topology_key]:
                    fwd_efficiency[topology_key][mode] = [0.0, 0.0]
                    fwd_events[topology_key][mode] = [0.0, 0.0]
                    num_samples[topology_key][mode] = [0, 0]

                events = (event_probs[EVENT_LLM] - 1) + event_probs[EVENT_SLM]

                # if len(fwd_events[topology_key][mode]) < 2:
                #     fwd_events[topology_key][mode].append(events)
                #     fwd_efficiency[topology_key][mode].append(float(path_lengths[topology_key]) / float(events))
                #     num_samples[topology_key][mode].append(1)
                # else:
                fwd_events[topology_key][mode][annc_sizes.index(annc_size)] += events
                # print("fwd_efficiency %s:%s:%d : %f" % (topology_key, mode, mode, (float(path_lengths[topology_key]) / float(events))))
                fwd_efficiency[topology_key][mode][annc_sizes.index(annc_size)] += (float(4) / float(events))
                num_samples[topology_key][mode][annc_sizes.index(annc_size)] += 1

                print("num_samples %s:%d:%d : %s" % (topology_key, mode, annc_size, str(num_samples[topology_key][mode])))

            if file_type == "path":

                # this will save the average delivery latency value
                avg_delivery_latency = 0.0

                delivery_deductable = [0.0, 0.0]
                path = get_path(os.path.join(test_dir, ("%s.test" % (topology_key))), file_label)
                df_aux = data[file_type][file_label].groupby(by = ["AS","STATUS"])["PROB"].sum()

                if 0 in df_aux[path[2]]:
                    delivery_deductable[0] = df_aux[path[2]][0]
                if 1 in df_aux[path[3]]:
                    delivery_deductable[1] = df_aux[path[3]][1]

                # keep track of the delivery probability
                delivery_prob = 0.0
                delivery_probs = data[file_type][file_label].groupby(by = ["STATUS"])["PROB"].sum()

                if topology_key not in avg_delivery_latencies:
                    avg_delivery_latencies[topology_key] = defaultdict()
                    avg_delivery_probs[topology_key] = defaultdict()

                if mode not in avg_delivery_latencies[topology_key]:
                    avg_delivery_latencies[topology_key][mode] = [0.0, 0.0]
                    avg_delivery_probs[topology_key][mode] = [[0.0, 0.0], [0.0, 0.0]]

                avg_delivery_probs[topology_key][mode][annc_sizes.index(annc_size)][0] += delivery_probs[0] - delivery_deductable[0]
                avg_delivery_probs[topology_key][mode][annc_sizes.index(annc_size)][1] += delivery_probs[1] - delivery_deductable[1]

                # order dataframe values by STATUS and LATENCY. the objective is 
                # to only look at the STATUS = 0 (correct delivery) and keep 
                # multiplying LATENCY * PROB until we reach a delivery prob of 1.0
                delivery_latency = data[file_type][file_label].sort(["STATUS","LATENCY"])
                
                for index, row in delivery_latency.iterrows():

                    print("%s.%d.%d : %d, %d, %.2f" % (topology_key, mode, annc_size, row['STATUS'], row['LATENCY'], row['PROB']))
                    if (delivery_prob + row['PROB']) >= 1.0:
                        avg_delivery_latency += ((1.0 - delivery_prob) * row['LATENCY'])
                        break

                    delivery_prob += row['PROB']
                    avg_delivery_latency += (row['PROB'] * row['LATENCY'])

                avg_delivery_latencies[topology_key][mode][annc_sizes.index(annc_size)] += avg_delivery_latency

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

    print("avg_delivery_latencies :")
    for t in avg_delivery_latencies:
        for b in avg_delivery_latencies[t]:
            print("%s.%s = %s" % (t, b, str(avg_delivery_latencies[t][b])))

    print("avg_delivery_probs :")
    for t in avg_delivery_probs:
        for b in avg_delivery_probs[t]:
            print("%s.%s = %s" % (t, b, str(avg_delivery_probs[t][b])))

    # print("path_lengths :")
    # for t in path_lengths:
    #     print("%s : %d" % (t, path_lengths[t]))

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
    bar_width = 0.20

    x_pos = defaultdict(list)
    xx_pos = defaultdict(list)

    show_legend = True
    i = 0
    for mode in modes:

        if mode != 1:
            show_legend = False

        for t, topology_key in enumerate(topology_keys):

            if show_legend:
                leg = topology_labels[topology_key]
            else:
                leg = None

            xx_pos[i] = (np.arange(1, (2 * len(modes)), step = 2) + (m * bar_width) + (bar_width / 2.0))
            i += 1
            ax1.bar(np.arange(1, (2 * len(modes)), step = 2) + (m * bar_width), np.array(avg_delivery_latencies[topology_key][mode]) / np.array(num_samples[topology_key][mode]), color = topology_colors[t], linewidth = 1.5, alpha = 0.55, width = bar_width, label = leg)
            m += 1.0

        x_pos[mode] = np.arange(1, (2 * len(modes)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    ax1.set_xlabel("Multiple match res. mode\nCache annc. size")
    ax1.set_ylabel("Avg. latency (# of hops)")
    # ax1.set_yscale('log')
    ax1.set_ylim(0, 8)
    ax1.set_yticks([0, 2, 4, 6, 8])

    # xticks = x_pos['256']
    xxticks = x_pos[1] + ((x_pos[2][0] - x_pos[1][0]) / 2.0)
    xticks = interleave_n(x_pos[1], xxticks, x_pos[2])
    # a convoluted way to set the x-axis labels?    
    xtick_labels = []
    for annc_size in annc_sizes:
        xtick_labels.append("%s" % ('RML'))
        xtick_labels.append("\n%d" % (annc_size))
        xtick_labels.append("%s" % ('Fallbk'))

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
        for mode in modes:
            for topology_key in topology_keys:
                avg_correct_delivery_probs_array.append(float(avg_delivery_probs[topology_key][mode][f][0] / num_samples[topology_key][mode][f]))
                avg_incorrect_delivery_probs_array.append(float(avg_delivery_probs[topology_key][mode][f][1] / num_samples[topology_key][mode][f]))

    # another convoluted way to create the x-axis for the fwd efficiency graph
    # FIXME: interleave_n() accepts a list of lists
    print(xx_pos)
    xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5])
    # xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2])

    ax2.plot(xx_pos, avg_correct_delivery_probs_array, linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = 'Corr. del.')
    ax2.plot(xx_pos, avg_incorrect_delivery_probs_array, linewidth = 1.5, color = 'black', linestyle = '--', markersize = 5, marker = 'v', label = 'Incorr. del.')
    ax2.axhspan(0, 2, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    ax2.set_xlim(xx_pos[0] - (3 * bar_width), xx_pos[-1] + (3 * bar_width))
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)

    ax2.set_ylim(-2, 6)
    # ax2.set_yscale('log')
    ax2.set_yticks([0, 1, 2])
    ax2.set_ylabel("Avg. # of deliveries")

    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("cdn-fallbacks.pdf", bbox_inches='tight', format = 'pdf')

def plot_fallbacks_no_tps(data, test_dir):

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
    bf_keys = ['192', '256']

    # modes : '1' relf, '2' fallback
    modes = [1, 2]

    # annc. size
    annc_sizes = [2, 5]

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
            if bf_key == '192':
                continue

            # collect the mode as int
            mode = int(file_label.split("-")[5])
            # if mode not in modes:
            #     modes.append(mode)
            annc_size = 0

            # if not ('fp2' in file_label or 'fp5' in file_label):
            #     continue

            if file_label.split("-")[7] != 'fp5':
                annc_size = 2
            else:
                annc_size = 5

            if file_type == "events":

                # gather the total probs by event
                event_probs = data[file_type][file_label].groupby(by = ["EVENT"])["PROB"].sum()

                # gather the probabilities for MLM and SLM events for 
                # forwarding efficiency (these events are mutually exclusive)
                if topology_key not in fwd_efficiency:
                    fwd_efficiency[topology_key] = defaultdict()
                    fwd_events[topology_key] = defaultdict()
                    num_samples[topology_key] = defaultdict()

                if mode not in fwd_efficiency[topology_key]:
                    fwd_efficiency[topology_key][mode] = [0.0, 0.0]
                    fwd_events[topology_key][mode] = [0.0, 0.0]
                    num_samples[topology_key][mode] = [0, 0]

                events = (event_probs[EVENT_LLM] - 1) + event_probs[EVENT_SLM]

                # if len(fwd_events[topology_key][mode]) < 2:
                #     fwd_events[topology_key][mode].append(events)
                #     fwd_efficiency[topology_key][mode].append(float(path_lengths[topology_key]) / float(events))
                #     num_samples[topology_key][mode].append(1)
                # else:
                fwd_events[topology_key][mode][annc_sizes.index(annc_size)] += events
                # print("fwd_efficiency %s:%s:%d : %f" % (topology_key, mode, mode, (float(path_lengths[topology_key]) / float(events))))
                fwd_efficiency[topology_key][mode][annc_sizes.index(annc_size)] += (float(4) / float(events))
                num_samples[topology_key][mode][annc_sizes.index(annc_size)] += 1

                print("num_samples %s:%s:%d : %s" % (topology_key, mode, mode, str(num_samples[topology_key][mode])))

            if file_type == "path":

                # this will save the average delivery latency value
                avg_delivery_latency = 0.0

                # hack warning : removes unnecessary correct and incorrect delivery probs, 
                # added to .tsv files by mistake
                delivery_deductable = [0.0, 0.0]
                path = get_path(os.path.join(test_dir, ("%s.test" % (topology_key))), file_label)
                df_aux = data[file_type][file_label].groupby(by = ["AS","STATUS"])["PROB"].sum()

                # if an incorrect delivery happens for a node 2 hops away from the 
                # request source, the last delivery (to the origin server) should be removed
                if 0 in df_aux[path[2]]:
                    delivery_deductable[0] = 1
                # this was supposed to compensate for a delivery on the last CDN in a path 
                # in-between request source and origin server
                if 1 in df_aux[path[2]]:
                    delivery_deductable[1] = df_aux[path[2]][1]

                delivery_probs = data[file_type][file_label].groupby(by = ["STATUS"])["PROB"].sum()

                if topology_key not in avg_delivery_latencies:
                    avg_delivery_latencies[topology_key] = defaultdict()
                    avg_delivery_probs[topology_key] = defaultdict()

                if mode not in avg_delivery_latencies[topology_key]:
                    avg_delivery_latencies[topology_key][mode] = [0.0, 0.0]
                    avg_delivery_probs[topology_key][mode] = [[0.0, 0.0], [0.0, 0.0]]

                avg_delivery_probs[topology_key][mode][annc_sizes.index(annc_size)][0] += delivery_probs[0] - delivery_deductable[0]
                avg_delivery_probs[topology_key][mode][annc_sizes.index(annc_size)][1] += delivery_probs[1] - delivery_deductable[1]

                # order dataframe values by STATUS and LATENCY. the objective is 
                # to only look at the STATUS = 0 (correct delivery) and keep 
                # multiplying LATENCY * PROB until we reach a delivery prob of 1.0.
                # FIXME : this seems wrong, and should be changed
                delivery_latency = data[file_type][file_label].sort(["STATUS","LATENCY"])

                # keep track of the delivery probability
                delivery_prob = 0.0
                for index, row in delivery_latency.iterrows():

                    print("%s.%s.%d : %d, %d, %.2f" % (topology_key, mode, mode, row['STATUS'], row['LATENCY'], row['PROB']))
                    if (delivery_prob + row['PROB']) >= 1.0:
                        avg_delivery_latency += ((1.0 - delivery_prob) * row['LATENCY'])
                        break

                    delivery_prob += row['PROB']
                    avg_delivery_latency += (row['PROB'] * row['LATENCY'])

                avg_delivery_latencies[topology_key][mode][annc_sizes.index(annc_size)] += avg_delivery_latency

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

    print("avg_delivery_latencies :")
    for t in avg_delivery_latencies:
        for b in avg_delivery_latencies[t]:
            print("%s.%s = %s" % (t, b, str(avg_delivery_latencies[t][b])))

    print("avg_delivery_probs :")
    for t in avg_delivery_probs:
        for b in avg_delivery_probs[t]:
            print("%s.%s = %s" % (t, b, str(avg_delivery_probs[t][b])))

    # print("path_lengths :")
    # for t in path_lengths:
    #     print("%s : %d" % (t, path_lengths[t]))

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
    bar_width = 0.20

    x_pos = defaultdict(list)
    xx_pos = defaultdict(list)

    show_legend = True
    i = 0
    for mode in modes:

        if mode != 1:
            show_legend = False

        for t, topology_key in enumerate(topology_keys):

            if show_legend:
                leg = topology_labels[topology_key]
            else:
                leg = None

            xx_pos[i] = (np.arange(1, (2 * len(modes)), step = 2) + (m * bar_width) + (bar_width / 2.0))
            i += 1
            ax1.bar(np.arange(1, (2 * len(modes)), step = 2) + (m * bar_width), np.array(avg_delivery_latencies[topology_key][mode]) / np.array(num_samples[topology_key][mode]), color = topology_colors[t], linewidth = 1.5, alpha = 0.55, width = bar_width, label = leg)
            m += 1.0

        x_pos[mode] = np.arange(1, (2 * len(modes)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    ax1.set_xlabel("Multiple match res. mode\nCache annc. size")
    ax1.set_ylabel("Avg. latency (# of hops)")
    # ax1.set_yscale('log')
    ax1.set_ylim(0, 14)
    ax1.set_yticks([0, 2, 4, 6, 8, 10])

    # xticks = x_pos['256']
    xxticks = x_pos[1] + ((x_pos[2][0] - x_pos[1][0]) / 2.0)
    xticks = interleave_n(x_pos[1], xxticks, x_pos[2])
    # a convoluted way to set the x-axis labels?    
    xtick_labels = []
    for annc_size in annc_sizes:
        xtick_labels.append("%s" % ('RML'))
        xtick_labels.append("\n%d" % (annc_size))
        xtick_labels.append("%s" % ('Fallbk'))

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
        for mode in modes:
            for topology_key in topology_keys:
                avg_correct_delivery_probs_array.append(float(avg_delivery_probs[topology_key][mode][f][0] / num_samples[topology_key][mode][f]))
                avg_incorrect_delivery_probs_array.append(float(avg_delivery_probs[topology_key][mode][f][1] / num_samples[topology_key][mode][f]))

    # another convoluted way to create the x-axis for the fwd efficiency graph
    # FIXME: interleave_n() accepts a list of lists
    print(xx_pos)
    xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5])
    # xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2])

    ax2.plot(xx_pos, avg_correct_delivery_probs_array, linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = 'Corr. del.')
    ax2.plot(xx_pos, avg_incorrect_delivery_probs_array, linewidth = 1.5, color = 'black', linestyle = '--', markersize = 5, marker = 'v', label = 'Incorr. del.')
    ax2.axhspan(0, 3, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    ax2.set_xlim(xx_pos[0] - (3 * bar_width), xx_pos[-1] + (3 * bar_width))
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)

    ax2.set_ylim(-1, 6)
    # ax2.set_yscale('log')
    ax2.set_yticks([0, 1, 2, 3])
    ax2.set_ylabel("Avg. # of deliveries")

    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("cdn-fallbacks-no-tps.pdf", bbox_inches='tight', format = 'pdf')