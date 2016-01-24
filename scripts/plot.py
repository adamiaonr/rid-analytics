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

AVERAGE_LATENCY="average_latency"
CACHE_LATENCY="cache_latency"
ORIGIN_LATENCY="origin_latency"

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 3:
        print "usage: python plot.py <input-data-filename>.csv <output-filename>.png"
        return

    # extract the .csv file argument
    data_file = sys.argv[1];

    if not os.path.isfile(data_file):
        print "ERROR: argument " + data_file + " doesn't exist. abort."
        return

    # extract the .csv file components into a list of tuples which will contain 
    # the values <latency, probability, decision_type>, converting strings to
    # floats
    f = open(data_file, 'rb')

    latency_data_list = []

    # we will draw 2 vertical lines for reference: latency to cache and origin 
    # server
    average_latency = 0.0
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
        elif splitted[0] == AVERAGE_LATENCY:
            average_latency = float(splitted[1])
        else:
            latency_data_list.append((float(splitted[0]), float(splitted[1]), str(splitted[2]).rstrip()))

    # get index of repeated latencies in latency_data_list
    x = [ seq[0] for seq in latency_data_list ]
    repeated = []

    for i in xrange(len(x) - 1):
        if x[i + 1] == x[i]:
            repeated.append(i);

    # condense items with the same latency
    c = collections.Counter()
    for k, v1, v2 in latency_data_list:
        c[k] += v1 

    # let's plot the damn thing then...
    fig = plt.figure()
    cdf_plot = fig.add_subplot(111)

    # sort c by key
    c = sorted(c.items())
    x = [ seq[0] for seq in c ]
    y = [ seq[1] for seq in c ]

    y = np.cumsum(y)

    # make sure everything's less than 1 (approx. errors)
    for i in xrange(len(y)):
        if y[i] > 1.0:
            y[i] = 1.0

    prev_x = min(x)
    prev_y = min(y)

    cdf_plot.grid(True)

    y_max = 1.0

    log_scale = 0

    if (log_scale):
        cdf_plot.set_yscale('log')
        y_min = math.pow(10, math.floor(math.log10(min(y))))
    else:
        y_min = 0.0

    # plot avg. latency line
    average_latency_line = cdf_plot.plot([average_latency, average_latency], [y_min, y_max], 'b--', linewidth=1.5)
    # plot cache latency line
    cache_latency_line = cdf_plot.plot([cache_latency, cache_latency], [y_min, y_max], 'g--', linewidth=1.5)
    # plot origin server latency
    origin_latency_line = cdf_plot.plot([origin_latency, origin_latency], [y_min, y_max], 'r--', linewidth=1.5)    

    # plot the discrete CDF by parts
    for i in xrange(len(x)):

        cdf_plot.plot([prev_x, x[i]], [prev_y, prev_y], 'b-', linewidth=1.5)

        # including the vertical lines between 
        cdf_plot.plot([x[i], x[i]], [prev_y, y[i]], 'b-', linewidth=.50)
        prev_x = x[i]
        prev_y = y[i]

    cdf_plot.plot(x[1:], y[:-1], 'wo', ms=8)
    cdf_plot.plot(x, y, 'ro', ms=10)

    cdf_plot.set_xlabel('Latency (hops)', fontsize=LABEL_FONT_SIZE)
    cdf_plot.set_ylabel('CDF', fontsize=LABEL_FONT_SIZE)
    cdf_plot.legend((average_latency_line[0], cache_latency_line[0], origin_latency_line[0]), ('average latency', 'cache latency', 'orig. server latency'), loc='upper left', fontsize=LEGEND_FONT_SIZE)
    
    # append the avg, cache and origin latencies to x
    x.append(cache_latency)
    x.append(average_latency)
    x.append(origin_latency)

    cdf_plot.set_xlim(0.0, custom_ceil(max(x)))
    cdf_plot.set_xticks(np.arange(0.0, custom_ceil(max(x)) + 1, 1.0))
    cdf_plot.set_ylim(y_min, y_max)

#    plt.show()
    plt.savefig(sys.argv[2], bbox_inches='tight')


if __name__ == "__main__":
    main()
 

