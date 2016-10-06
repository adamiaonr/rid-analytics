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

event_labels = {0: 'No matches', 1: 'Multiple links', 2: 'Local link', 3: 'Single link'}

def extract_data(data_dir):

    data = defaultdict(OrderedDict)

    for file_name in sorted(glob.glob(os.path.join(data_dir, '*.tsv'))):

        file_type = file_name.split(".")[0].split("/")[-1]
        file_label = file_name.split(".")[1]

        print("type = %s, label = %s" % (file_type, file_label))

        data[file_type][file_label] = pd.read_csv(file_name, sep = "\t")
        data[file_type][file_label] = data[file_type][file_label].convert_objects(convert_numeric = True)

    return data

def plot_base(data):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    event_avgs = defaultdict(list)
    path_values = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No matches', 'Single link', 'Multiple links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe', 'white']
    # labels and styles for lines (path probs)
    path_labels = ['Correct deliv.', 'Wrong deliv.']
    path_colors = ['black', 'black']
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

    ax1.set_title("Avg. nr. of events & path outcome probs.\n(|R|=15, 192 bit BF)")
    ax1.set_xlabel("Fwd. entry size", fontsize = 12)
    ax1.set_ylabel("Avg. nr. of events", fontsize = 12)
    ax1.set_xticks(np.arange(1, (2 * len(entry_lengths)), step = 2), entry_lengths)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(np.arange(1, (2 * len(entry_lengths)), step = 2), np.array(path_values[l]), linewidth = 1.5, color = path_colors[l], marker = path_markers[l], label = path_labels[l])

    ax2.set_ylim(0.0, 2.5)
    ax2.set_ylabel("Path outcome prob.")
    ax2.set_xticks(np.arange(1, (2 * len(entry_lengths)), step = 2))
    ax2.set_xticklabels(entry_lengths)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("base.pdf", bbox_inches='tight', format = 'pdf')

def plot_bf_sizes(data):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    event_avgs = defaultdict(list)
    path_values = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No matches', 'Single link', 'Multiple links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe', 'white']
    # labels and styles for lines (path probs)
    path_labels = ['Correct deliv.', 'Wrong deliv.']
    path_colors = ['black', 'black']
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

    ax1.set_title("Avg. nr. of events & path outcome probs.\n(|R|=15, |F|=1)")
    ax1.set_xlabel("BF sizes (in bit)", fontsize = 12)
    ax1.set_ylabel("Avg. nr. of events", fontsize = 12)
    ax1.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2), bf_sizes)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(np.arange(1, (2 * len(bf_sizes)), step = 2), np.array(path_values[l]), linewidth = 1.5, color = path_colors[l], marker = path_markers[l], label = path_labels[l])

    ax2.set_ylim(0.0, 2.5)
    ax2.set_ylabel("Path outcome prob.")
    ax2.set_xticks(np.arange(1, (2 * len(bf_sizes)), step = 2))
    ax2.set_xticklabels(bf_sizes)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("bf-sizes.pdf", bbox_inches='tight', format = 'pdf')

def plot_req_sizes(data):

    # the point is to plot the avg. nrs. of events as bars, and the prob. of 
    # path outcomes as a line (because these are diff. things)
    event_avgs = defaultdict(list)
    path_values = defaultdict(list)

    # labels and colors for bars (events)
    event_labels = ['No matches', 'Single link', 'Multiple links', 'Local link']
    event_colors = ['#000000', '#708090', '#bebebe', 'white']
    # labels and styles for lines (path probs)
    path_labels = ['Correct deliv.', 'Wrong deliv.']
    path_colors = ['black', 'black']
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

    ax1.set_title("Avg. nr. of events & path outcome probs.\n(|F|=1, 256 bit BF)")
    ax1.set_xlabel("Req. size", fontsize = 12)
    ax1.set_ylabel("Avg. nr. of events", fontsize = 12)
    ax1.set_xticks(np.arange(1, (2 * len(req_sizes)), step = 2), req_sizes)
    ax1.legend(fontsize=12, ncol=1, loc='upper left')

    # line plot on top of bar chart
    ax2 = ax1.twinx()

    for l in np.arange(len(path_labels)):
        ax2.plot(np.arange(1, (2 * len(req_sizes)), step = 2), np.array(path_values[l]), linewidth = 1.5, color = path_colors[l], marker = path_markers[l], label = path_labels[l])

    ax2.set_ylim(0.0, 2.5)
    ax2.set_ylabel("Path outcome prob.")
    ax2.set_xticks(np.arange(1, (2 * len(req_sizes)), step = 2))
    ax2.set_xticklabels(req_sizes)
    ax2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
    ax2.set_yticklabels(['0.0', '0.5', '1.0', '', '', ''])
    ax2.legend(fontsize=12, ncol=1, loc='upper right')

    plt.savefig("req-sizes.pdf", bbox_inches='tight', format = 'pdf')

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
        plot_base(data)
    elif args.case == 'bf-sizes':
        plot_bf_sizes(data)
    elif args.case == 'req-sizes':
        plot_req_sizes(data)
    else:
        sys.stderr.write("""%s: [ERROR] please supply a valid case\n""" % sys.argv[0]) 
        parser.print_help()
        sys.exit(1)