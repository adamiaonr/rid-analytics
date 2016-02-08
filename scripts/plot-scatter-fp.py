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

    if len(sys.argv) < 3:
        print "usage: python plot-scatter-fp.py <input-file-base-name> <output-filename>.png"
        return

    # extract the .csv file argument
    base_name = sys.argv[1]

    fig = plt.figure()
    subplot_code = 220

    for i in xrange(4):

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

        # ax.set_ylim(0.0, max(d[max(d, key=lambda i: d[i])]) + 1)
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
        
        origin_latency = 7.0
        origin_latency_line = ax.axhline(y=origin_latency, c='r', ls='--', linewidth=2.0)
        # plot cache latency line
        cache_latency = 5.0
        cache_latency_line = ax.axhline(y=cache_latency, c='g', ls='--', linewidth=2.0)

        if (i == 0):
            a = plt.plot([], [], color='pink', linewidth=10)
            b = plt.plot([], [], color='lightgreen', linewidth=10)

            #ax.legend((a[0], b[0], origin_latency_line, cache_latency_line), ('Feedback', 'Fallback', 'To origin', 'To cache'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

        ax.set_xlabel("False positive rate (@tier " + str(i + 1) + ")", fontsize=LABEL_FONT_SIZE)
        ax.set_ylabel('Avg. Latency (nr. of hops)', fontsize=LABEL_FONT_SIZE)
        ax.tick_params(axis='x', labelsize=10)

    plt.savefig(sys.argv[2], bbox_inches='tight')

#         feedback_latencies = []

#         for line in f.readlines():

#             splitted = line.split(",")        
#             feedback_latencies.append((float(splitted[0]), float(splitted[1])))

#         _x = [ seq[0] for seq in feedback_latencies ]
#         _y = [ seq[1] for seq in feedback_latencies ]

#         _y_min = math.floor(min(_y)) - 1
#         _y_max = math.ceil(max(_y)) + 1

#         #a = ax.scatter(_x, _y, c='c', s=100, alpha=0.50)
#         a = ax.boxplot(_y, positions = [1, 2], widths = 0.6)

#         filename = base_name + "." + str(i + 1) + ".fallback.csv"
#         print filename

#         if not os.path.isfile(filename):
#             print "ERROR: argument " + filename + " doesn't exist. abort."
#             return

#         ax = fig.add_subplot(subplot_code)

#         f = open(filename, 'rb')

#         fallback_latencies = []

#         for line in f.readlines():

#             splitted = line.split(",")        
#             fallback_latencies.append((float(splitted[0]), float(splitted[1])))

#         x = [ seq[0] for seq in fallback_latencies ]
#         y = [ seq[1] for seq in fallback_latencies ]

#         y_min = math.floor(min(y)) - 1
#         y_max = math.ceil(max(y)) + 1

#         b = ax.scatter(x, y, c='m', s=50, alpha=0.50)    

#         # plot avg. latency line
#         origin_latency = 7.0
#         origin_latency_line = ax.axhline(y=origin_latency, c='r', ls='--', linewidth=2.0)
#         # plot cache latency line
#         cache_latency = 1.0
#         cache_latency_line = ax.axhline(y=cache_latency, c='g', ls='--', linewidth=2.0)

#         ax.set_xlabel("False positive rate (@tier " + str(i + 1) + ")", fontsize=LABEL_FONT_SIZE)
#         ax.set_ylabel('Avg. latency (nr. of hops)', fontsize=LABEL_FONT_SIZE)

#         if i == 0:
#             ax.legend((a, b, origin_latency_line, cache_latency_line), ('Feedback', 'Fallback', 'Hops to origin', 'Hops to cache'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

#         ax.grid(True)
#         fig.tight_layout()
        
#         ax.set_xscale('log')
#         stackd.set_xticks([0.0000000001, 0.00000001, 0.0000001, 0.0001, 0.1, 1.0])
#         ax.set_xlim(0.0000000001, 1.0)
#         ax.set_ylim(min([y_min, _y_min]), max([y_max, _y_max]))

# #    plt.show()
#     plt.savefig(sys.argv[2], bbox_inches='tight')

if __name__ == "__main__":
    main()