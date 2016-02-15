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

LABEL_FONT_SIZE=12
LEGEND_FONT_SIZE=12

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 6:
        print "usage: python plot-box-fp.py <nr of tiers> <cache lat.> <origin lat.> <input-file-base-name> <output-filename>.png"
        return

    # extract the .csv file argument
    nr_tiers = int(sys.argv[1])
    cache_latency = int(sys.argv[2])
    origin_latency = int(sys.argv[3])
    base_name = sys.argv[4]

    fig = plt.figure()
    subplot_code = (2 * 100) + (2 * 10)

    for i in xrange(nr_tiers):

        filename = base_name + "." + str(i + 1) + ".feedback.csv"

        if not os.path.isfile(filename):
            print "ERROR: argument " + filename + " doesn't exist. abort."
            return

        f = open(filename, 'rb')

        _temp_data = []

        for line in f.readlines():
            splitted = line.split(",")        
            _temp_data.append((float(splitted[0]), float(splitted[1])))

        data = defaultdict(list)
        for k, v in _temp_data:
            data[k].append(v)

        d = collections.OrderedDict(sorted(data.items()))

        subplot_code += 1
        ax = fig.add_subplot(subplot_code)

        xx = []
        for key in d.keys():
            if key < 0.5:
                xx.append(str("1E" + str(int(np.log10(key)))));
            else:
                xx.append(str(key))

        print xx

        print d.keys()
        x = np.arange(1.0, len(d.keys()) + 1)
        print (x - 0.5)

        bp = ax.boxplot(d.values(), positions = ((2.0 * x) - 0.25) - 1, widths=0.50, patch_artist=True);
        ax.set_xticks((2.0 * x) - 1)

        colors = ['pink']
        for patch in bp['boxes']:
            patch.set_facecolor(colors[0])

        ax.set_ylim((0.75 * float(cache_latency)), max(max(d[max(d, key=lambda i: d[i])]) + 1, 1.25 * float(origin_latency)))
        # ax.set_yticks(np.arange(0.0, max(d[max(d, key=lambda i: d[i])]) + 2));

        ####

        filename = base_name + "." + str(i + 1) + ".fallback.csv"

        if not os.path.isfile(filename):
            print "ERROR: argument " + filename + " doesn't exist. abort."
            return

        f = open(filename, 'rb')

        _temp_data = []

        for line in f.readlines():
            splitted = line.split(",")        
            _temp_data.append((float(splitted[0]), float(splitted[1])))

        data = defaultdict(list)
        for k, v in _temp_data:
            data[k].append(v)

        d = collections.OrderedDict(sorted(data.items()))

        print d.keys()
        x = np.arange(1.0, len(d.keys()) + 1)
        print (x + 0.5)

        bp = ax.boxplot(d.values(), positions = (((2.0 * x) + 0.25)) - 1, widths=0.50, patch_artist=True);
        ax.set_xticks((2.0 * x) - 1)
        ax.set_xticklabels(xx)


        colors = ['lightgreen']
        for patch in bp['boxes']:
            patch.set_facecolor(colors[0])

        ax.set_xlim(min(((2.0 * x) - 0.25) - 1) - 1, max((((2.0 * x) + 0.25)) - 1) + 1)
        
        origin_latency_line = ax.axhline(y=origin_latency, c='r', ls='--', linewidth=2.0)
        cache_latency_line = ax.axhline(y=cache_latency, c='g', ls='--', linewidth=2.0)

        ax.set_xlabel("False positive rate (@tier " + str(i + 1) + ")", fontsize=LABEL_FONT_SIZE)
        ax.set_ylabel('Avg. Latency (nr. of hops)', fontsize=LABEL_FONT_SIZE)
        ax.tick_params(axis='x', labelsize=10)

        if (i == (nr_tiers - 1)):

            subplot_code += 1
            ax = fig.add_subplot(subplot_code)
            
            a = plt.plot([], [], color='pink', linewidth=10)
            b = plt.plot([], [], color='lightgreen', linewidth=10)
            ax.legend((a[0], b[0], origin_latency_line, cache_latency_line), ('Feedback', 'Fallback', 'To origin', 'To cache'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

    plt.savefig(sys.argv[5], bbox_inches='tight')

if __name__ == "__main__":
    main()