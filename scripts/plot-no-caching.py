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

# in honor of "stack'd", the pittsburgh restaurant next to my shadyside 
# place

LABEL_FONT_SIZE=12
LEGEND_FONT_SIZE=10

BAR_WIDTH=0.35

# outcome type strings
outcomes = ['Correct', 'Incorrect', 'Fallback', 'Relay', 'TP', 'FP', 'Unknown']
colors = ['green', 'red', 'orange', 'black']

# sub dirs with results for different request sizes
sub_dirs = ['r5', 'r10']
request_sizes = [5, 10]

# distributions of P(|F|_i)
f_distributions = []

for i in xrange(len(sub_dirs)):
    f_distributions.append([])

    for j in xrange(2):
        f_distributions[i].append([])

f_distributions[0][0].append([0.0004,0.0458,1.3624,17.141,81.45])
f_distributions[0][1].append([81.45,17.147,1.3552,0.0465,0.001])
f_distributions[1][0].append([0.0001,0.0001,0.0001,0.0001,0.0025,0.0593,0.7756,6.289,29.88,62.994])
f_distributions[1][1].append([62.99,29.876,6.2915,0.7736,0.0646,0.004,0.0001,0.0001,0.0001,0.0001])

# distributions of P(|F\R|)
f_r_distributions = []

for i in xrange(len(sub_dirs)):
    f_r_distributions.append([])

    for j in xrange(2):
        f_r_distributions[i].append([])

f_r_distributions[0][0].append([0.0004,0.0458,1.3624,17.141,81.45])
f_r_distributions[0][1].append([81.45,17.147,1.3552,0.0465,0.001])
f_r_distributions[1][0].append([0.0001,0.0001,0.0001,0.0001,0.0025,0.0593,0.7756,6.289,29.88,62.994])
f_r_distributions[1][1].append([62.99,29.876,6.2915,0.7736,0.0646,0.004,0.0001,0.0001,0.0001,0.0001])

print f_distributions

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-outcomes.py <bf-size> <input-file-dir> <output-file-dir> <output-file-name>"
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

        _for_real_dir = input_file_dir + "/" + sub_dir

        for file_name in os.listdir(_for_real_dir):

            if file_name.endswith(bf_size + ".csv"):

                _file = open(_for_real_dir + "/" + file_name, 'rb')

                for line in _file.readlines():
                    line_splitted = line.split(",") 

                    try:
                        outcome_index = int(int(line_splitted[4]))
                        request_size = ((int(file_name.split("-")[2]) / 5) - 1)
                        origin_distance = int(file_name.split("-")[1])
                        _data[request_size][outcome_index].append((origin_distance, float(line_splitted[2])))

                        if (float(line_splitted[2]) > 0.0):

                            _latency_to_origin = float(line_splitted[3])
                            _latency[request_size].append([origin_distance, _latency_to_origin * float(line_splitted[2])])
                    
                    except IndexError:

                        print "line = " + line
                        print "out. index = " + str(outcome_index) + " out. str = " + outcomes[outcome_index]
                        print "request size = " + str(request_size) + ", " + file_name.split("-")[2] + ", " + str((int(file_name.split("-")[2]) / 5) - 1)
                        print "origin distance = " + str(origin_distance)
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

    # subplot_code = (2 * 100) + (2 * 10)
    # distributions of P(|F|_i)
    # subplot_code += 1
    # f_distributions_plot = fig.add_subplot(subplot_code)

    # markers for this graph
    markers = ['--o', '-d', ':o', ':x']

    for i in xrange(len(sub_dirs)):

        f_distributions_plot = plt.subplot2grid((4,4), (i,0), colspan=2)

        # a grid makes it easier to read
        f_distributions_plot.grid(True)

        # title goes here
        f_distributions_plot.set_title("|R| = 5", fontsize=LABEL_FONT_SIZE)

        for j in xrange(2):

            f_distributions_plot.plot((np.arange(request_sizes[i]) + 1), f_distributions[i][j][0], markers[j], linewidth=2, color='black')

        # labels
        if (i > 0):
            f_distributions_plot.set_xlabel("Prefix length |F|", fontsize=LEGEND_FONT_SIZE)
        f_distributions_plot.set_ylabel('%', fontsize=LABEL_FONT_SIZE)

        f_distributions_plot.set_ylim(0.0, 100.0)

        f_distributions_plot.set_xticks(np.arange(1, request_sizes[i] + 1, 1))
        x_labels = f_distributions_plot.get_xticklabels()
        plt.setp(x_labels, rotation=0, fontsize=LEGEND_FONT_SIZE)

        y_labels = f_distributions_plot.get_yticklabels()
        plt.setp(y_labels, rotation=0, fontsize=LEGEND_FONT_SIZE)

        _mlocal = plt.plot([], [], '--o', color='black', linewidth=2)
        _m_non_local = plt.plot([], [], '-d', color='black', linewidth=2)

        f_distributions_plot.legend((_mlocal[0], _m_non_local[0]), ('local iface', 'non-local iface'), loc='upper center', fontsize=LEGEND_FONT_SIZE)

    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
    plt.savefig(sys.argv[3] + "/" + sys.argv[4] + "-dist" + ".pdf", format='pdf', bbox_inches='tight')

    # # 2 legends
    # _mlocal = plt.plot([], [], '-o', color='black', linewidth=2)
    # _m_non_local = plt.plot([], [], '-x', color='black', linewidth=2)

    # legend1 = f_distributions_plot.legend((_r5[0], _r10[0]), ('|R| = 5', '|R| = 10'), loc=1, fontsize=LEGEND_FONT_SIZE)

    # f_distributions_plot.add_artist(legend1)

    # markers for this graph
    markers = ['-^', '-o', '-d', '-v']
    markers_size = [7.5, 10, 7.5, 7.5]

    for i in xrange(len(sub_dirs)):

        # subplot_code += 1
        # outcome_probs = fig.add_subplot(subplot_code)

        fig = plt.figure()
        outcome_probs = plt.subplot2grid((6,6), (0,0), rowspan=2, colspan=2)

        # outcome_probs.set_title("Final state probabilities for |R| = " + str(request_sizes[i]), fontsize=LABEL_FONT_SIZE)
        outcome_probs.grid(True)

        for j in xrange(len(outcomes)):

            if (j > 3):
                break

            outcome_probs.plot(x[i][j], y[i][j], markers[j], markersize=markers_size[j], linewidth=2, color=colors[j])

            outcome_probs.set_xlabel("Origin dist. (# of hops)", fontsize=LABEL_FONT_SIZE)
            outcome_probs.set_ylabel('Outcome prob. [0.0, 1.0]', fontsize=LABEL_FONT_SIZE)

            outcome_probs.set_yscale('log')

        _cd = plt.plot([], [], markers[0], markersize=5, color=colors[0], linewidth=2)
        _id = plt.plot([], [], markers[1], markersize=5, color=colors[1], linewidth=2)
        _fd = plt.plot([], [], markers[2], markersize=5, color=colors[2], linewidth=2)
        _rd = plt.plot([], [], markers[3], markersize=5, color=colors[3], linewidth=2)

        outcome_probs.legend((_cd[0], _id[0], _fd[0], _rd[0]), (outcomes[0], outcomes[1], outcomes[2], outcomes[3]), loc='center left', fontsize=LEGEND_FONT_SIZE)

        plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
        plt.savefig(sys.argv[3] + "/" + sys.argv[4] + "-outcomes-" + sub_dirs[i] + ".pdf", format='pdf', bbox_inches='tight')

    # just an experiment with stack'd graphs...
    # for i in xrange(len(sub_dirs)):

    #     if (i < 1):
    #         continue

    #     # here, we cumulatively sum the columns of the data[] array : the idea is 
    #     # to sum the elements below row i to the row i, which results in an offset 
    #     # for the data in row i (this needs a better explanation)
    #     y_stack = np.cumsum(y[i], axis = 0)
    #     y_stack = (y_stack / max(y_stack[3])) * 100.0

    #     # subplot_code += 1
    #     # stackd = fig.add_subplot(subplot_code)
    #     stackd = plt.subplot2grid((4,4), (2,0), rowspan=2, colspan=2)

    #     stackd.set_title("% of final states for |R| = " + str(request_sizes[i]), fontsize=LABEL_FONT_SIZE)
    #     stackd.grid(True)

    #     _x_min = 1000
    #     _x_max = 0

    #     for j in xrange(4):

    #         if (j > 3):
    #             break

    #         if (max(x[i][j]) > _x_max):
    #             _x_max = max(x[i][j])

    #         if (min(x[i][j]) < _x_min):
    #             _x_min = min(x[i][j])

    #         if (j == 0):
    #             stackd.fill_between(x[i][j], 0, y_stack[1,:], facecolor=colors[j], alpha=.7)
    #         else:
    #             stackd.fill_between(x[i][j], y_stack[j - 1,:], y_stack[j,:], facecolor=colors[j], alpha=.7)

    #     stackd.set_xlabel('Origin distance (# of hops)', fontsize=LABEL_FONT_SIZE)
    #     stackd.set_ylabel('Final state %', fontsize=LABEL_FONT_SIZE)
        
    #     _cd = plt.plot([], [], '-', color=colors[0], linewidth=5)
    #     _id = plt.plot([], [], '-', color=colors[1], linewidth=5)
    #     _fd = plt.plot([], [], '-', color=colors[2], linewidth=5)
    #     _rd = plt.plot([], [], '-', color=colors[3], linewidth=5)

    #     stackd.legend((_cd[0], _id[0], _fd[0], _rd[0]), (outcomes[0], outcomes[1], outcomes[2], outcomes[3]), loc='upper right', fontsize=LEGEND_FONT_SIZE)

    #     stackd.grid(True)

    #     # x-axis handling by parts:
    #     stackd.set_xlim(_x_min, _x_max)

    #     # y axis is easy...
    #     stackd.set_ylim(0.0, 100.0)

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

    avg_latencies.plot(np.arange(1, 10, 1), np.arange(1, 10, 1), '--', linewidth=2, color='black');

    avg_latencies.set_xlim(0, 10)
    avg_latencies.set_xticks(np.arange(1, 10, 1))

    avg_latencies.set_ylim(0, math.ceil(y_max + 1))
    avg_latencies.set_yticks(np.arange(1, math.ceil(y_max + 1), 2))

    avg_latencies.set_xlabel('Origin dist. (# of hops)', fontsize=LABEL_FONT_SIZE)
    avg_latencies.set_ylabel('Avg. latency (# of hops)', fontsize=LABEL_FONT_SIZE)

    _r5 = plt.plot([], [], '-', color=_colors[0], linewidth=5)
    _r10 = plt.plot([], [], '-', color=_colors[1], linewidth=5)
    _ori = plt.plot([], [], '--', color='green', linewidth=2)

    avg_latencies.legend((_r5[0], _r10[0]), ('|R| = 5', '|R| = 10'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
    plt.savefig(sys.argv[3] + "/" + sys.argv[4] + "-latencies" + ".pdf", format='pdf', bbox_inches='tight')

if __name__ == "__main__":
    main()
