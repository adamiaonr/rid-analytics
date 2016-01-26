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
LEGEND_FONT_SIZE=10

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 3:
        print "usage: python plot-scatter.py <input-file-base-name> <output-filename>.png"
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

        subplot_code += 1
        ax = fig.add_subplot(subplot_code)

        f = open(filename, 'rb')

        feedback_latencies = []

        for line in f.readlines():

            splitted = line.split(",")        
            feedback_latencies.append((float(splitted[0]), float(splitted[1])))

        _x = [ seq[0] for seq in feedback_latencies ]
        _y = [ seq[1] for seq in feedback_latencies ]

        _y_min = math.floor(min(_y)) - 1
        _y_max = math.ceil(max(_y)) + 1

        a = ax.scatter(_x, _y, c='c', s=100, alpha=0.50)

        filename = base_name + "." + str(i + 1) + ".fallback.csv"
        print filename

        if not os.path.isfile(filename):
            print "ERROR: argument " + filename + " doesn't exist. abort."
            return

        ax = fig.add_subplot(subplot_code)

        f = open(filename, 'rb')

        fallback_latencies = []

        for line in f.readlines():

            splitted = line.split(",")        
            fallback_latencies.append((float(splitted[0]), float(splitted[1])))

        x = [ seq[0] for seq in fallback_latencies ]
        y = [ seq[1] for seq in fallback_latencies ]

        y_min = math.floor(min(y)) - 1
        y_max = math.ceil(max(y)) + 1

        b = ax.scatter(x, y, c='m', s=50, alpha=0.50)    

        # plot avg. latency line
        origin_latency = 10.0
        origin_latency_line = ax.axhline(y=origin_latency, c='r', ls='--', linewidth=2.0)
        # plot cache latency line
        cache_latency = 7.0
        cache_latency_line = ax.axhline(y=cache_latency, c='g', ls='--', linewidth=2.0)

        ax.set_xlabel("false positive rate (level " + str(i + 1) + ")", fontsize=LABEL_FONT_SIZE)
        ax.set_ylabel('avg. latency (hops)', fontsize=LABEL_FONT_SIZE)

        if i == 0:
            ax.legend((a, b, origin_latency_line, cache_latency_line), ('feedback', 'fallback', 'to origin', 'to cache'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

        ax.grid(True)
        fig.tight_layout()
        
        ax.set_xscale('log')
        ax.set_xlim(0.0000000001, 1.0)
        ax.set_ylim(min([y_min, _y_min]), max([y_max, _y_max]))

#    plt.show()
    plt.savefig(sys.argv[2], bbox_inches='tight')

if __name__ == "__main__":
    main()