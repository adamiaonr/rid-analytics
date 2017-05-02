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

def get_bar_labels(path_lengths):

    length_labels = defaultdict()

    for length in path_lengths:

        # labels should be in the format:
        # "<Length 1>: <min-length1>-<max-length1> hop"
        # "<Length 2>: <min-length2>"" (in case min == max)

        # determine the max and min lengths
        max_length = 0
        min_length = 256

        for outdegree in path_lengths[length]:
            if path_lengths[length][outdegree] > max_length:
                max_length = path_lengths[length][outdegree]
            if path_lengths[length][outdegree] < min_length:
                min_length = path_lengths[length][outdegree]

        # capitalize the length label (e.g. 'short' becomes 'Short')
        label_prefix = length.capitalize()
        # if the top label, add a 'path' suffix to it, to make it more 
        # descriptive
        if length == 'short':
            label_prefix = label_prefix + " path"

        # distinguish between min == max, min != max
        if max_length == min_length:
            length_labels[length] = ("%s: %d" % (label_prefix, min_length))
        else:
            length_labels[length] = ("%s: %d-%d" % (label_prefix, min_length, max_length))

        # add a 'hop' suffix to the top label, to make it more 
        # descriptive
        if length == 'short':
            length_labels[length] += " hop"

    return length_labels

def plot_efficiency(data, test_file):

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

    # labels and colors for bars (path lengths)
    length_colors = ['#000000', '#708090', '#bebebe']
    length_keys = ['short', 'median', 'long']
    outdegree_keys = ['low', 'median']

    # labels and styles for lines (path probs)
    delivery_labels = ['Corr. del.', 'Wrong del.']
    delivery_colors = ['black', 'black', 'black', 'black']
    delivery_styles = ['-', '--', '-.', ':']
    delivery_markrs = ['o', 'v', '^', 'o']
    delivery_markrs_size = [10, 10, 10, 8]

    # x-axis labels (|f| values)
    entry_lengths = []

    # populate the fwd_* dicts w/ data obtained from .tsv files
    got_path_info = False
    for file_type in data:
        for file_label in data[file_type]:

            if not got_path_info:
                path_lengths, avg_outdegrees = get_path_info(test_file, file_label)
                got_path_info = True

            entry_length = int(file_label.split("-")[3])
            if entry_length not in entry_lengths:
                # print("scanning |f| = %d" % (entry_length))
                entry_lengths.append(entry_length)

            # collect path lengths and avg. outdegrees
            length_key = file_label.split("-")[7]
            outdegree_key = file_label.split("-")[8]

            if outdegree_key not in outdegree_keys:
                continue

            if file_type == "events":

                # gather the total probs by event
                event_probs = data[file_type][file_label].groupby(by = ["EVENT"])["PROB"].sum()

                # gather the probabilities for MLM and SLM events for 
                # forwarding efficiency (these events are mutually exclusive)
                if length_key not in fwd_efficiency:
                    fwd_efficiency[length_key] = defaultdict()
                    fwd_events[length_key] = defaultdict()

                if outdegree_key not in fwd_efficiency[length_key]:
                    fwd_efficiency[length_key][outdegree_key] = []
                    fwd_events[length_key][outdegree_key] = []

                fwd_events[length_key][outdegree_key].append(event_probs[EVENT_MLM] + event_probs[EVENT_SLM])
                fwd_efficiency[length_key][outdegree_key].append(float(path_lengths[length_key][outdegree_key]) / float(event_probs[EVENT_MLM] + event_probs[EVENT_SLM]))

    # print("fwd_efficiency : %s" % (fwd_efficiency))
    # print("fwd_events :")
    # for l in fwd_events:
    #     for o in fwd_events[l]:
    #         print("%s.%s = %s" % (l, o, str(fwd_events[l][o])))

    # create labels for the bars
    length_labels = get_bar_labels(path_lengths)

    #matplotlib.style.use('ggplot')
    fig = plt.figure(figsize=(5, 3.5))
    ax1 = fig.add_subplot(111)
    ax1.xaxis.grid(False)
    ax1.yaxis.grid(True)

    # each |f| slot will have 3 groups of bars, w/ 
    # 3 bars each. there will be a total of 
    # 27 bars in the graph
    bar_group_size = len(length_keys)
    bar_group_num = 2

    # assumes the inter bar group space is half a bar. also, for n groups of bars 
    # we have n - 1 inter bar group spaces
    m = -(float(bar_group_num * bar_group_size) / 1.0) - ((bar_group_num - 1) / 2.0)
    bar_width = 0.20

    x_pos = defaultdict(list)
    xx_pos = defaultdict(list)

    show_legend = True
    i = 0
    for outdegree in outdegree_keys:

        if outdegree != 'low':
            show_legend = False

        for l, length in enumerate(length_keys):

            if show_legend:
                leg = length_labels[length]
            else:
                leg = None

            if outdegree not in fwd_events[length]:
                fwd_events[length][outdegree] = [0, 0, 0]

            xx_pos[i] = (np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width) + (bar_width / 2.0))
            i += 1
            ax1.bar(np.arange(1, (2 * len(entry_lengths)), step = 2) + (m * bar_width), np.array(fwd_events[length][outdegree]), color = length_colors[l], linewidth = 1.5, alpha = 0.55, width = bar_width, label = leg)
            m += 1.0

        x_pos[outdegree] = np.arange(1, (2 * len(entry_lengths)), step = 2) + ((m - 1.5) * bar_width)
        m += 1.0

    ax1.set_xlabel("Avg. path outdegree\nFwd. entry size")
    ax1.set_ylabel("Avg. # of used links")
    # ax1.set_yscale('log')
    ax1.set_ylim(0, 40)
    ax1.set_yticks([0, 10, 20, 30])

    # a convoluted way to set the x-axis labels?
    xpto = defaultdict(list)
    for length in path_lengths:
        for outdegree in path_lengths[length]:
            if outdegree not in avg_outdegrees[length]:
                continue
            xpto[outdegree].append(avg_outdegrees[length][outdegree])

    xxticks = x_pos['low'] + ((x_pos['median'][0] - x_pos['low'][0]) / 2.0)
    xticks = interleave_n(x_pos['low'], xxticks, x_pos['median'])
    # a convoluted way to set the x-axis labels?
    xtick_labels = []
    for entry_length in entry_lengths:
        xtick_labels.append("%.1f" % (np.mean(xpto['low'])))
        xtick_labels.append("\n%d" % (entry_length))
        xtick_labels.append("%.1f" % (np.mean(xpto['median'])))

    # for e in entry_lengths[1:]:
    #     xtick_labels.append('')
    #     xtick_labels.append('\n%d' % (e))
    #     xtick_labels.append('')

    ax1.set_xticks(xticks)
    ax1.set_xticklabels(xtick_labels, fontsize = 14)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # fwd_efficiency as a line plot on top of bar chart
    ax2 = ax1.twinx()
    ax2.yaxis.grid(False)
    # a convoluted way to transform a dictionary into a list, 
    # ready to use in plot(). the outermost guide of the cycle 
    # is the entry size |f|
    fwd_efficiency_array = []
    for f, el in enumerate(entry_lengths):
        for outdegree in outdegree_keys:
            for length in ['short', 'median', 'long']:

                if outdegree not in fwd_efficiency[length]:
                    fwd_efficiency[length][outdegree] = [-1.0, -1.0, -1.0]

                fwd_efficiency_array.append(int(fwd_efficiency[length][outdegree][f] * 100.0))

    # another convoluted way to create the x-axis for the fwd efficiency graph
    # FIXME: interleave_n() accepts a list of lists
    print(xx_pos)
    xx_pos = interleave_n(xx_pos[0], xx_pos[1], xx_pos[2], xx_pos[3], xx_pos[4], xx_pos[5])

    # remove the indexes from xx_pos and fwd_efficiency_array which have 
    # no values
    indices_to_remove = [i for i, x in enumerate(fwd_efficiency_array) if x == -100]
    for i in sorted(indices_to_remove, reverse = True):
        del fwd_efficiency_array[i]
        del xx_pos[i]

    ax2.plot(xx_pos, fwd_efficiency_array, linewidth = 1.5, color = 'black', linestyle = '-', markersize = 5, marker = 'o', label = 'Fwd. effic.')
    # ax2.axhspan(0, 100, linewidth = 0.0, facecolor = '#bebebe', alpha=0.20)

    ax2.set_xlim(xx_pos[0] - (2 * bar_width), xx_pos[-1] + (2 * bar_width))
    ax2.set_xticks(xticks)
    ax2.set_xticklabels(xtick_labels)

    ax2.set_ylim(-50, 150)
    ax2.set_yticks([0, 50, 100])
    ax2.set_ylabel("Fwd. efficiency (%)")

    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("base-case-efficiency.pdf", bbox_inches='tight', format = 'pdf') 