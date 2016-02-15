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
LEGEND_FONT_SIZE=12

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-stackd.py <nr. of tiers> <input-file-base-name> <output-filename>.png"
        return

    # extract the .csv file argument
    base_name = sys.argv[2]
    nr_tiers = int(sys.argv[1])

    outcomes = ['cc', 'ic', 'rlyd', 'drpd']

    data = [[], [], [], []]

    fig = plt.figure()
    subplot_code = (2 * 100) + (2 * 10)

    for n in xrange(nr_tiers):
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

        stackd.fill_between(x, 0, y_stack[0,:], facecolor='lightgreen', alpha=.7)
        stackd.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor='pink', alpha=.7)
        stackd.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor='cyan', alpha=.7)
        stackd.fill_between(x, y_stack[2,:], y_stack[3,:], facecolor='gray', alpha=.7)

        stackd.set_xlabel("False positive rate (@tier " + str(n + 1) + ")", fontsize=LABEL_FONT_SIZE)
        stackd.set_ylabel('Outcome (%)', fontsize=LABEL_FONT_SIZE)
        
        a = plt.plot([], [], color='lightgreen', linewidth=10)
        b = plt.plot([], [], color='pink', linewidth=10)
        c = plt.plot([], [], color='cyan', linewidth=10)
        d = plt.plot([], [], color='gray', linewidth=10)

        stackd.grid(True)
        fig.tight_layout()

        stackd.set_xscale('log')
        stackd.set_xticks([0.000000001, 0.0000001, 0.00001, 0.001, 0.1, 1.0])
        stackd.set_xlim(0.000000001, 1.0)
        stackd.set_ylim(0.0, 100.0)

        if (n == (nr_tiers - 1)):
            subplot_code += 1
            stackd = fig.add_subplot(subplot_code)
            stackd.legend((a[0], b[0], c[0], d[0]), ('Correct dest.', 'Wrong dest.', 'Relayed', 'Dropped'), loc='lower right', fontsize=LEGEND_FONT_SIZE)

#    plt.show()
    plt.savefig(sys.argv[3], bbox_inches='tight')

if __name__ == "__main__":
    main()