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
outcomes_short = ['CD', 'Incorrect delivery', 'Fallback delivery', 'Fallback relay', 'TP', 'FP', 'Unknown']
colors = ['green', 'red', 'orange', 'black']

exp_factors = [1, 2, 5, 10]

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-prefix-management.py <bf-size> <input-file-dir> <output-file-dir> <output-file-name>"
        return

    start_time = time.time()

    bf_size = sys.argv[1]
    input_file_dir = sys.argv[2]

    _data = []

    for i in xrange(len(outcomes)):
        _data.append([])

        for j in xrange(len(exp_factors)):
            _data[i].append([])

    # read each .csv file in input_file_dir
    _real_dir = input_file_dir + "/r10"
    _limit = 7

    for file_name in os.listdir(_real_dir):

        if file_name.endswith(".csv"):

            _file = open(_real_dir + "/" + file_name, 'rb')

            for line in _file.readlines():
                line_splitted = line.split(",") 

                try:
                    outcome_index = int(int(line_splitted[4]))
                    request_size = ((int(file_name.split("-")[2]) / 5) - 1)
                    origin_distance = int(file_name.split("-")[1])
                    exp_factor = int(line_splitted[6])
                    f_min = int(line_splitted[5])

                    if (f_min > _limit):
                        break

                    _data[outcome_index][exp_factors.index(exp_factor)].append((f_min, float(line_splitted[2])))
                
                except IndexError:

                    print "line = " + line
                    print "out. index = " + str(outcome_index) + " out. str = " + outcomes[outcome_index]
                    print "request size = " + str(request_size) + ", " + file_name.split("-")[2] + ", " + str((int(file_name.split("-")[2]) / 5) - 1)
                    print "origin distance = " + str(origin_distance)
                    print "exp_factor = " + str(exp_factor)
                    print "f_min = " + str(f_min)
                    print "prob. = " + str(float(line_splitted[2]))

                    return

    elapsed_time = time.time() - start_time
    print "[READ FILES IN " + str(elapsed_time) + " sec]"

    x = []
    y = []

    for i in xrange(len(outcomes)):
        x.append([])
        y.append([])

        if (i > 3):
            break

        for j in xrange(len(exp_factors)):

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

    subplot_code = (2 * 100) + (2 * 10)

    # markers for this graph
    markers = ['-o', '-v', '--x', ':d']

    for i in xrange(len(outcomes)):

        if (i > 3):
            break

    fig = plt.figure()

    avg_latencies = plt.subplot2grid((4,4), (0,0), rowspan=2, colspan=2)

        # a grid makes it easier to read
        prefix_man_plot.grid(True)

        # title goes here
        if (i < 1):
            prefix_man_plot.set_title("Correct del. (160 bit, |R| = 10)", fontsize=LABEL_FONT_SIZE)
        else:
            prefix_man_plot.set_title(outcomes_short[i], fontsize=LABEL_FONT_SIZE)

        for j in xrange(len(exp_factors)):

            prefix_man_plot.plot(x[i][j], y[i][j], markers[j], linewidth=2, color=colors[i])

        # labels
        prefix_man_plot.set_xlabel("Min. entry length |F|min", fontsize=LEGEND_FONT_SIZE)
        prefix_man_plot.set_ylabel("Outcome probability [0.0, 1.0]", fontsize=LEGEND_FONT_SIZE)

        _exp_a = plt.plot([], [], '-o', color=colors[i], linewidth=2)
        _exp_b = plt.plot([], [], '-v', color=colors[i], linewidth=2)
        _exp_c = plt.plot([], [], '--x', color=colors[i], linewidth=2)
        _exp_d = plt.plot([], [], ':d', color=colors[i], linewidth=2)

        _y_ticks = [1E-15, 1E-10, 1E-5, 1E0]

        if (i == 0):
            prefix_man_plot.set_yscale('log')
            prefix_man_plot.set_ylim(1E-15, 1.0)
            prefix_man_plot.set_yticks(_y_ticks)
            prefix_man_plot.legend((_exp_a[0], _exp_b[0], _exp_c[0], _exp_d[0]), (str(exp_factors[0]), str(exp_factors[1]), str(exp_factors[2]), str(exp_factors[3])), loc='lower left', fontsize=LEGEND_FONT_SIZE)

        if (i == 1):
            prefix_man_plot.set_yscale('log')
            prefix_man_plot.set_ylim(1E-15, 1.0)
            prefix_man_plot.set_yticks(_y_ticks)
            prefix_man_plot.legend((_exp_a[0], _exp_b[0], _exp_c[0], _exp_d[0]), (str(exp_factors[0]), str(exp_factors[1]), str(exp_factors[2]), str(exp_factors[3])), loc='lower left', fontsize=LEGEND_FONT_SIZE)

        if (i == 2):
            prefix_man_plot.set_ylim(0.6, 1.0)
            prefix_man_plot.legend((_exp_a[0], _exp_b[0], _exp_c[0], _exp_d[0]), (str(exp_factors[0]), str(exp_factors[1]), str(exp_factors[2]), str(exp_factors[3])), loc='lower right', fontsize=LEGEND_FONT_SIZE)

        if (i == 3):
            prefix_man_plot.set_yscale('log')
            prefix_man_plot.set_ylim(1E-15, 1.0)
            prefix_man_plot.set_yticks(_y_ticks)
            prefix_man_plot.legend((_exp_a[0], _exp_b[0], _exp_c[0], _exp_d[0]), (str(exp_factors[0]), str(exp_factors[1]), str(exp_factors[2]), str(exp_factors[3])), loc='lower left', fontsize=LEGEND_FONT_SIZE)


        prefix_man_plot.set_xticks(np.arange(1, _limit + 1, 1))


    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=0.4)
    plt.savefig(sys.argv[3] + "/" + sys.argv[4] + ".pdf", format='pdf', bbox_inches='tight')

if __name__ == "__main__":
    main()
