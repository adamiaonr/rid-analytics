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

# in honor of "Stack'd", the pittsburgh restaurant next to my shadyside 
# place

LABEL_FONT_SIZE=10
LEGEND_FONT_SIZE=10

outcomes = ['cc', 'ic', 'rlyd', 'drpd']
penalty_types = ['Feedback', 'Fallback']

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-stackd.py <input-file-dir> <output-file-dir>"
        return

    # extract the .csv file argument
    input_file_dir = sys.argv[1]

    data = [[], [], [], []]

    fig = plt.figure()
    subplot_code = (2 * 100) + (2 * 10)

    for penalty in penalty_types:
        for i in xrange(len(outcomes)):

            filename = input_file_dir + "/outcomes." + penalty.lower() + "." + outcomes[i] + ".csv"

            if not os.path.isfile(filename):
                print "ERROR: " + filename + " doesn't exist. abort."
                return

            # open the file for a particular outcome & penalty combination
            f = open(filename, 'rb')

            # extract data from it
            _data = []

            for line in f.readlines():
                line_splitted = line.split(",")        
                _data.append((str(line_splitted[0]).upper(), float(line_splitted[3])))

            # a Counter allows us to 'reduce' all entries with similar 
            # <fp>.<alpha> pairs 
            _counter = collections.Counter()
            for k, v in _data:
                _counter[k] += v

            _counter = sorted(_counter.items())

            x = [ seq[0] for seq in _counter ]
            data[i] = [ seq[1] for seq in _counter ]

        print "X AXIS VALUES: " + str(x)

        # accumulate the values from diff. outcomes, for the same <fp>.<alpha>
        # pair 
        y_stack = np.cumsum(data, axis = 0)
        # make them represented from 0.0 to 100.0 %
        y_stack = (y_stack / max(y_stack[len(outcomes) - 1])) * 100.0
        print "Y AXIS VALUES: "
        print y_stack

        subplot_code += 1
        stackd = fig.add_subplot(subplot_code)

        _x = np.arange(len(x))

        stackd.fill_between(_x, 0, y_stack[0,:], facecolor='lightgreen', alpha=.7)
        stackd.fill_between(_x, y_stack[0,:], y_stack[1,:], facecolor='pink', alpha=.7)
        stackd.fill_between(_x, y_stack[1,:], y_stack[2,:], facecolor='cyan', alpha=.7)
        stackd.fill_between(_x, y_stack[2,:], y_stack[3,:], facecolor='gray', alpha=.7)

        stackd.set_xlabel("<FP Rate>.<Alpha> pairs", fontsize=LABEL_FONT_SIZE)
        stackd.set_ylabel('Outcome occurrence %', fontsize=LABEL_FONT_SIZE)
        
        a = plt.plot([], [], color='lightgreen', linewidth=10)
        b = plt.plot([], [], color='pink', linewidth=10)
        c = plt.plot([], [], color='cyan', linewidth=10)
        d = plt.plot([], [], color='gray', linewidth=10)

        stackd.grid(True)
        fig.tight_layout()

        # x-axis handling by parts:
        stackd.set_xticks(_x)
        stackd.set_xlim(min(_x), max(_x))
        stackd.set_xticklabels(x)   
        x_labels = stackd.get_xticklabels()
        plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

        # y axis is easy...
        stackd.set_ylim(0.0, 100.0)

        # title for subplot
        stackd.set_title(penalty)

        # legend
        stackd.legend((a[0], b[0], c[0], d[0]), ('Correct', 'Incorrect', penalty, 'Dropped'), loc='lower right', fontsize=LEGEND_FONT_SIZE)

    plt.savefig(sys.argv[2] + "/stackd.cache." + input_file_dir[-1:] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()