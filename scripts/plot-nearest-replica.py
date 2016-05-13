from datetime import date
from datetime import datetime
import matplotlib.pyplot as plt
import matplotlib.mlab as mlab
from matplotlib.pyplot import figure, show
from matplotlib import *
import matplotlib
import time
import sys
import os
import numpy as np
import collections
import math

LABEL_FONT_SIZE=12
LEGEND_FONT_SIZE=10

# width for bar char (latencies)
BAR_WIDTH=0.35

# outcome type strings
outcomes = ['Correct', 'Incorrect', 'Fallback', 'Relay', 'TP', 'FP', 'Unknown']
colors = ['green', 'red', 'orange', 'black']

# sub dirs with results for different request sizes
sub_dirs = ['r5', 'r10']
request_sizes = [5, 10]

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-nearest-replica.py <bf-size> <input-file-dir> <output-file-dir> <output-file-name>"
        return

    start_time = time.time()

    bf_size = sys.argv[1]
    input_file_dir = sys.argv[2]

    _data = []
    _latency = []

    # we will directly fill a data array structured as 
    # _data[<outcome index>][(<scenario code>, <probability>)]
    for i in xrange(len(sub_dirs)):
        _data.append([])
        _latency.append([])

        for j in xrange(len(outcomes)):
            _data[i].append([])

    # read each .csv file in input_file_dir
    for sub_dir in sub_dirs:

        real_dir = input_file_dir + "/" + sub_dir

        for file_name in os.listdir(real_dir):

            if file_name.endswith(bf_size + ".csv"):

                data_file = open(real_dir + "/" + file_name, 'rb')

                for line in data_file.readlines():
                    line_splitted = line.split(",") 

                    try:
                        outcome_index = int(int(line_splitted[4]))
                        request_size = ((int(file_name.split("-")[3]) / 5) - 1)
                        cache_distance = int(file_name.split("-")[2])
                        _data[request_size][outcome_index].append((cache_distance, float(line_splitted[2])))

                        if (float(line_splitted[2]) > 0.0):

                            latency_to_origin = float(line_splitted[3])
                            _latency[request_size].append([cache_distance, latency_to_origin * float(line_splitted[2])])
                    
                    except IndexError:

                        print "line = " + line
                        print "outcome index = " + str(outcome_index) + " outcome str = " + outcomes[outcome_index]
                        print "request size = " + str(request_size) + ", " + file_name.split("-")[3] + ", " + str((int(file_name.split("-")[3]) / 5) - 1)
                        print "cache distance = " + str(cache_distance)
                        print "prob. = " + str(float(line_splitted[2]))

                        return

    elapsed_time = time.time() - start_time
    print "[READ FILES IN " + str(elapsed_time) + " sec]"

    x = []
    y = []

    for i in xrange(len(sub_dirs)):
        x.append([])
        y.append([])

        print sub_dirs[i]

        for j in xrange(len(outcomes)):

            if (j > 3):
                break

            x[i].append([])
            y[i].append([])

            # this sums all probabilities of the same <outcome, scn code> pair 
            # into a single bucket data[outcome][scn code]
            _counter = collections.Counter()
            for key, value in _data[i][j]:
                _counter[key] += value

            # sort [scn codes] alphabetically
            _counter = sorted(_counter.items())

            x[i][j] = [ seq[0] for seq in _counter ]
            y[i][j] = [ seq[1] for seq in _counter ]

            print x[i][j]
            print y[i][j]

    # now lets generate a 2 x 2 plot, i.e. with 4 subplots, one for each data series
    fig = plt.figure()

    subplot_code = (3 * 100) + (3 * 10)

    # markers for this graph
    markers = ['-^', '-o', '-d', '-v']
    markers_size = [7.5, 10, 7.5, 7.5]

    for i in xrange(len(sub_dirs)):

        fig = plt.figure()
        outcome_probs = plt.subplot2grid((6,6), (0,0), rowspan=2, colspan=2)

        # outcome_probs.set_title("|R| = " + str(request_sizes[i]) + ", BF size " + bf_size + "bit", fontsize=LABEL_FONT_SIZE)
        outcome_probs.grid(True)

        for j in xrange(len(outcomes)):

            if (j > 3):
                break

            outcome_probs.plot(x[i][j], y[i][j], markers[j], markersize=markers_size[j], linewidth=2, color=colors[j])

            outcome_probs.set_xlabel('Cache dist. (# of hops)', fontsize=LABEL_FONT_SIZE)

            outcome_probs.set_ylabel('Outcome prob. [0.0, 1.0]', fontsize=LABEL_FONT_SIZE)
            outcome_probs.set_yscale('log')

            if (i == 0):
                outcome_probs.set_ylim(1E-12, 1E0)
                outcome_probs.set_yticks([1E-12, 1E-10, 1E-8, 1E-6, 1E-4, 1E-2, 1E0])

        _cd = plt.plot([], [], markers[0], markersize=5, color=colors[0], linewidth=2)
        _id = plt.plot([], [], markers[1], markersize=5, color=colors[1], linewidth=2)
        _fd = plt.plot([], [], markers[2], markersize=5, color=colors[2], linewidth=2)
        _rd = plt.plot([], [], markers[3], markersize=5, color=colors[3], linewidth=2)

        outcome_probs.legend((_cd[0], _id[0], _fd[0], _rd[0]), (outcomes[0], outcomes[1], outcomes[2], outcomes[3]), loc='lower center', fontsize=LEGEND_FONT_SIZE)

        plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
        plt.savefig(sys.argv[3] + "/" + sys.argv[4] + "-outcomes-" + sub_dirs[i] + ".pdf", format='pdf', bbox_inches='tight')

    x = []
    y = []

    bar_shift = -BAR_WIDTH

    fig = plt.figure()
    avg_latencies = plt.subplot2grid((6,6), (0,0), rowspan=2, colspan=2)
    avg_latencies.grid(True)
    # avg_latencies.set_title("Avg. latencies", fontsize=LABEL_FONT_SIZE)

    _colors = ['black', 'white']

    y_max = 0

    for i in xrange(len(sub_dirs)):

        print sub_dirs[i]

        x.append([])
        y.append([])

        # this sums all probabilities of the same <outcome, scn code> pair 
        # into a single bucket data[outcome][scn code]
        _counter = collections.Counter()
        for key, value in _latency[i]:
            _counter[key] += value

        # sort [scn codes] alphabetically
        _counter = sorted(_counter.items())

        x[i] = [ seq[0] for seq in _counter ]
        y[i] = [ seq[1] for seq in _counter ]

        for j in xrange(len(x[i])):
            x[i][j] += bar_shift

        print x[i]
        print y[i]

        if (y_max < max(y[i])):
            y_max = max(y[i])

        avg_latencies.bar(x[i], y[i], BAR_WIDTH, color=_colors[i])

        bar_shift = 0

    avg_latencies.plot(np.arange(1, 10, 1), np.arange(1, 10, 1), '--', linewidth=2, color='green');

    avg_latencies.set_xlim(0, 10)
    avg_latencies.set_xticks(np.arange(1, 10, 1))

    avg_latencies.set_ylim(0, math.ceil(y_max + 1))
    avg_latencies.set_yticks(np.arange(1, math.ceil(y_max + 1) + 1, 2))

    avg_latencies.set_xlabel('Cache dist. (# of hops)', fontsize=LABEL_FONT_SIZE)
    avg_latencies.set_ylabel('Avg. latency (# of hops)', fontsize=LABEL_FONT_SIZE)

    _r5 = plt.plot([], [], '-', color=_colors[0], linewidth=5)
    _r10 = plt.plot([], [], '-', color=_colors[1], linewidth=5)
#    _ori = plt.plot([], [], '--', color='green', linewidth=2)

    avg_latencies.legend((_r5[0], _r10[0]), ('|R| = 5', '|R| = 10'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
    plt.savefig(sys.argv[3] + "/" + sys.argv[4] + "-latencies" + ".pdf", format='pdf', bbox_inches='tight')

if __name__ == "__main__":
    main()
