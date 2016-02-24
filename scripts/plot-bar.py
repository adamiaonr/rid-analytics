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
from collections import defaultdict
import math
from pylab import plot, show, savefig, xlim, figure, \
                hold, ylim, legend, boxplot, setp, axes

LABEL_FONT_SIZE=10
LEGEND_FONT_SIZE=10

BAR_WIDTH=0.35

outcomes = ['cc', 'ic', 'rlyd', 'drpd']
penalty_types = ['Feedback', 'Fallback']


def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-bar.py <input-file-dir> <output-file-dir>"
        return

    # extract the .csv file argument
    input_file_dir = sys.argv[1]

    # misc variables
    colors = ['pink', 'lightgreen']
    colors_index = 0
    bar_shift = -BAR_WIDTH
    y_max = 0

    fig = plt.figure()
    subplot_code = (2 * 100) + (2 * 10)
    subplot_code += 1
    bar = fig.add_subplot(subplot_code)

    for penalty in penalty_types:

        filename = input_file_dir + "/latencies." + penalty.lower() + ".csv"

        if not os.path.isfile(filename):
            print "ERROR: " + filename + " doesn't exist. abort."
            return

        # open the file for a particular outcome & penalty combination
        f = open(filename, 'rb')

        # extract data from it
        _data = []

        for line in f.readlines():
            line_splitted = line.split(",")        
            _data.append((str(line_splitted[0]).upper(), float(line_splitted[2])))

        data = defaultdict(list)
        for k, v in _data:
            data[k].append(v)

        # we now have an ordered dict of key-value pairs
        d = collections.OrderedDict(sorted(data.items()))

        x = d.keys()
        _x = np.arange(0.0, len(x))

        # extract the y values as a list (format of OrderedDict not suitable 
        # for bar plot)
        y = []
        for v in d.values():
            y.append(v[0])

        if max(y) > y_max:
            y_max = max(y)

        bar.grid(True)

        avg_latencies = bar.bar(_x + bar_shift, y, BAR_WIDTH, color=colors[colors_index])

        # update bar shift (change signal)
        bar_shift = 0
        colors_index += 1

        bar.set_xlabel("<FP Rate>.<Alpha> pairs", fontsize=LABEL_FONT_SIZE)
        bar.set_ylabel('Avg. Latency (nr. of hops)', fontsize=LABEL_FONT_SIZE)

        # x axis is the most complicated
        bar.set_xlim(min(_x) - 1, max(_x) + 1)
        bar.set_xticks(_x)
        bar.set_xticklabels(x)
        x_labels = bar.get_xticklabels()
        plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

    a = plt.plot([], [], color='pink', linewidth=10)
    b = plt.plot([], [], color='lightgreen', linewidth=10)

    # plot cache & origin latency lines
    c = bar.axhline(y=((int(input_file_dir[-1:]) * 2) + 1), c='g', ls='--', linewidth=2.0)
    d = bar.axhline(y=((3 * 2) + 1), c='r', ls='--', linewidth=2.0)

    # set y axis limits after drawing the rest of the graph
    bar.set_ylim(0, math.ceil(1.1 * y_max))
    bar.set_yticks(np.arange(0.0, math.ceil(1.1 * y_max) + 1))

    # legend
    bar.legend((a[0], b[0], c, d), ('Feedback', 'Fallback', 'To cache', 'To origin'), loc='lower right', fontsize=LEGEND_FONT_SIZE)

    # save figure as bar.cache.<dist-of-cache>.png in the graph dir folder
    plt.savefig(sys.argv[2] + "/bar.cache." + input_file_dir[-1:] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()