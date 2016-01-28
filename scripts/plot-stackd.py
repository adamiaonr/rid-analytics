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

# in the honor of "Stack'd", the pittsburgh restaurant next to my shadyside 
# place...

LABEL_FONT_SIZE=16
LEGEND_FONT_SIZE=10

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-stackd.py <input-file-base-name> <output-filename>.png"
        return

    # extract the .csv file argument
    base_name = sys.argv[1]

    outcomes = ['wds', 'cc', 'ds']

    data = [[], [], []]

    fig = plt.figure()
    subplot_code = 220

    for n in xrange(4):
        for i in xrange(len(outcomes)):

            filename = base_name + "." + str(n + 1) + ".outcomes." + outcomes[i] + ".csv"

            if not os.path.isfile(filename):
                print "ERROR: argument " + filename + " doesn't exist. abort."
                return

            f = open(filename, 'rb')

            _temp_data = []

            for line in f.readlines():
                splitted = line.split(",")        
                _temp_data.append((float(splitted[0]), float(splitted[2])))

            _temp_counter = collections.Counter()
            for k, v in _temp_data:
                _temp_counter[k] += v

            _temp_counter = sorted(_temp_counter.items())

            x = [ seq[0] for seq in _temp_counter ]
            data[i] = [ seq[1] for seq in _temp_counter ]

        print x
        y_stack = np.cumsum(data, axis = 0)
        print y_stack
        y_stack = (y_stack / max(y_stack[len(outcomes) - 1])) * 100.0

        subplot_code += 1
        stackd = fig.add_subplot(subplot_code)

        stackd.fill_between(x, 0, y_stack[0,:], facecolor="#CC6666", alpha=.7)
        stackd.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor="#1DACD6", alpha=.7)
        stackd.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor="#6E5160")

        stackd.set_xlabel("false positive rate (hop #" + str(n + 1) + ")", fontsize=LABEL_FONT_SIZE)
        stackd.set_ylabel('outcome %', fontsize=LABEL_FONT_SIZE)
        
        a = plt.plot([], [], color="#CC6666", linewidth=10)
        b = plt.plot([], [], color="#1DACD6", linewidth=10)
        c = plt.plot([], [], color="#6E5160", linewidth=10)

        stackd.legend((a[0], b[0], c[0]), ('wrong dest. > server', 'correct cache', 'dropped > server'), loc='upper left', fontsize=LEGEND_FONT_SIZE)

        stackd.grid(True)
        fig.tight_layout()

        stackd.set_xscale('log')
        stackd.set_xlim(0.000000001, 1.0)
        stackd.set_ylim(0.0, 100.0)

#    plt.show()
    plt.savefig(sys.argv[2], bbox_inches='tight')

if __name__ == "__main__":
    main()