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

LABEL_FONT_SIZE=16
LEGEND_FONT_SIZE=13

CACHE_LATENCY="cache_latency"
ORIGIN_LATENCY="origin_latency"

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 3:
        print "usage: python plotbar.py <input-data-filename>.csv <output-filename>.png"
        return

    # extract the .csv file argument
    data_file = sys.argv[1]

    if not os.path.isfile(data_file):
        print "ERROR: argument " + data_file + " doesn't exist. abort."
        return

    # extract the .csv file components into a list of tuples which will contain 
    # the values <latency, probability, decision_type>, converting strings to
    # floats
    f = open(data_file, 'rb')

    avg_latency_data_list = []

    # we will draw 2 vertical lines for reference: latency to cache and origin 
    # server
    cache_latency = 0.0
    origin_latency = 0.0

    for line in f.readlines():
        
        splitted = line.split(",")

        # extract the avg, cache and origin server latencies, if the specifiers 
        # appear
        if  splitted[0] == CACHE_LATENCY:
            cache_latency = float(splitted[1])
        elif splitted[0] == ORIGIN_LATENCY:
            origin_latency = float(splitted[1])
        else:
            # the average latencies for each of the scenarios
            avg_latency_data_list.append((str(splitted[0]), float(splitted[1])))

    # let's plot the damn thing then...
    fig = plt.figure()
    bar_plot = fig.add_subplot(111)

    tick_labels = [ seq[0] for seq in avg_latency_data_list ]
    #legend_labels = [ seq[2] for seq in avg_latency_data_list ]
    y = [ seq[1] for seq in avg_latency_data_list ]

    bar_plot.grid(True)

    y_max = custom_ceil(max(y)) + 1.0
    y_min = 0.0

    # plot cache latency line
    cache_latency_line = bar_plot.axhline(y=cache_latency, c='g', ls='--', linewidth=2.0)
    # plot origin server latency
    origin_latency_line = bar_plot.axhline(y=origin_latency, c='r', ls='--', linewidth=2.0)    

    ind = np.arange(len(tick_labels))
    width = 0.35
    avg_latencies = bar_plot.bar(ind, y, width, color='c')

    bar_plot.set_xlabel('Scenarios', fontsize=LABEL_FONT_SIZE)
    bar_plot.set_ylabel('Latency (hops)', fontsize=LABEL_FONT_SIZE)
    bar_plot.set_ylim(y_min, y_max)

    bar_plot.set_xticks(ind + (width / 2.0))
    bar_plot.set_xticklabels(tick_labels)
#    bar_plot.set_yticks(np.arange(y_min, y_max + 1.0, 1.0))

    bar_plot.legend((avg_latencies[0], cache_latency_line, origin_latency_line), ('avg. latency', 'cache latency', 'orig. latency'), loc='upper right', fontsize=LEGEND_FONT_SIZE)

#    plt.show()
    plt.savefig(sys.argv[2], bbox_inches='tight')

if __name__ == "__main__":
    main()
 

