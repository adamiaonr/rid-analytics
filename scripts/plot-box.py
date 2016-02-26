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
BAR_WIDTH=-0.50

penalty_types = ['feedback', 'fallback']
colors = ['pink', 'lightgreen']

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    start_time = time.time()

    if len(sys.argv) < 2:
        print "usage: python plot-box.py <x-param> <nr-of-tiers> <input-file-dir> <output-file-dir>"
        return

    # extract the .csv file argument
    x_param_str = sys.argv[1]

    if x_param_str == 'fprob':
        x_param = 0
        x_descr = "FP Rate at Tier "
    elif x_param_str == 'alpha':
        x_param = 2
        x_descr = "Alpha at Tier "
    else:
        print "ERROR: unknown param: " + x_param_str
        return

    tier_depth = int(sys.argv[2])
    input_file_dir = sys.argv[3]

    # calculate cache and origin latency
    cache_latency = (float(input_file_dir[-1:]) * 2.0) + 1.0
    origin_latency = float(((tier_depth * 2.0) + 1.0))

    y_max = 0

    filename = input_file_dir + "/latencies.csv"

    if not os.path.isfile(filename):
        print "ERROR: " + filename + " doesn't exist. abort."
        return

    # prepare the arrays which will be used as placeholders for data
    _data = []

    for i in xrange(tier_depth):
        _data.append([])

        for j in xrange(len(penalty_types)):
            _data[i].append([])

    # open the file for a particular outcome & penalty combination
    f = open(filename, 'rb')

    for line in f.readlines():
        
        line_splitted = line.split(",")

        try:
            _data[int(line_splitted[0 + x_param])][penalty_types.index(str(line_splitted[4]))].append((float(line_splitted[1 + x_param]), float(line_splitted[5])))
        
        except IndexError:

            print line
            print int(line_splitted[0 + x_param])
            print penalty_types.index(str(line_splitted[4]))

            return

    elapsed_time = time.time() - start_time

    print "[READ FILE IN " + str(elapsed_time) + " sec]"

    fig = plt.figure()
    subplot_code = (2 * 100) + (2 * 10)

    BAR_WIDTH=-0.50

    for i in xrange(tier_depth):

        subplot_code += 1
        box = fig.add_subplot(subplot_code)

        for j in xrange(len(penalty_types)):

            data = defaultdict(list)

            for k, v in _data[i][j]:
                data[k].append(v)

            # we now have an ordered dict of key-value pairs
            d = collections.OrderedDict(sorted(data.items()))

            # nicer way of representing powers of 10
            _x = []
            for key in d.keys():
                if key < 0.5:
                    _x.append(str("1E" + str(int(np.log10(key)))));
                else:
                    _x.append(str(key))
            print _x
            print d.keys()

            x = np.arange(1.0, len(d.keys()) + 1)
            print (x + BAR_WIDTH)

            avg_latencies = box.boxplot(d.values(), positions = ((2.0 * x) + (BAR_WIDTH / 2.0)) - 1, widths=BAR_WIDTH, patch_artist=True);

            # invert the signal of BAR_WIDTH
            BAR_WIDTH = -BAR_WIDTH
            
            for patch in avg_latencies['boxes']:
                patch.set_facecolor(colors[j])

            if max(d[max(d, key=lambda i: d[i])]) + 1.0 > y_max:
                y_max = max(d[max(d, key=lambda i: d[i])]) + 1.0

            print y_max

        origin_latency_line = box.axhline(y=origin_latency, c='r', ls='--', linewidth=2.0)
        cache_latency_line = box.axhline(y=cache_latency, c='g', ls='--', linewidth=2.0)

        box.grid(True)

        box.set_xlabel(x_descr + str(i), fontsize=LABEL_FONT_SIZE)
        box.set_ylabel('Avg. Latency (nr. of hops)', fontsize=LABEL_FONT_SIZE)
        box.tick_params(axis='x', labelsize=10)

        box.set_xticklabels(_x)
        box.set_xticks((2.0 * x) - 1)
        box.set_xlim(min(((2.0 * x) - (BAR_WIDTH / 2.0)) - 1.0) - 1.0, max((((2.0 * x) + (BAR_WIDTH / 2.0))) - 1.0) + 1.0)

        box.set_ylim((0.75 * cache_latency), max(y_max, 1.25 * origin_latency))        
        # box.set_yticks(np.arange(0.0, max(d[max(d, key=lambda i: d[i])]) + 2));

        a = plt.plot([], [], color='pink', linewidth=LEGEND_FONT_SIZE)
        b = plt.plot([], [], color='lightgreen', linewidth=LEGEND_FONT_SIZE)

    subplot_code += 1
    box = fig.add_subplot(subplot_code)
    box.legend((a[0], b[0], origin_latency_line, cache_latency_line), ('Feedback', 'Fallback', 'To origin', 'To cache'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

    # save figure
    plt.savefig(sys.argv[4] + "/box." + x_param_str + ".cache." + input_file_dir[-1:] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()