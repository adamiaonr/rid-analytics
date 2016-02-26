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

outcomes = ['correct dest.', 'wrong dest. > orig. server', 'fp detect. > orig. server', 'dropped']
penalty_types = ['feedback', 'fallback']

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    start_time = time.time()

    if len(sys.argv) < 2:
        print "usage: python plot-stackd.py <x-param> <nr-of-tiers> <input-file-dir> <output-file-dir>"
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

    filename = input_file_dir + "/outcomes.csv"

    if not os.path.isfile(filename):
        print "ERROR: " + filename + " doesn't exist. abort."
        return

    # prepare the arrays which will be used as placeholders for data
    _data = []

    for i in xrange(len(penalty_types)):
        _data.append([])

        for j in xrange(tier_depth):
            _data[i].append([])

            for k in xrange(len(outcomes)):
                _data[i][j].append([])

    # open the file for a particular outcome & penalty combination
    f = open(filename, 'rb')

    for line in f.readlines():
        
        line_splitted = line.split(",")

        try:
            _data[penalty_types.index(str(line_splitted[4]))][int(line_splitted[0 + x_param])][outcomes.index(str(line_splitted[5]))].append((float(line_splitted[1 + x_param]), float(line_splitted[6])))
        except IndexError:

            print line
            print penalty_types.index(str(line_splitted[4]))
            print int(line_splitted[0 + x_param])
            print outcomes.index(str(line_splitted[5]))

            return

    elapsed_time = time.time() - start_time

    print "[READ FILE IN " + str(elapsed_time) + " sec]"

    # done with the file reading. draw stack'd graphs for each penalty case.
    for i in xrange(len(penalty_types)):

        fig = plt.figure()
        subplot_code = (2 * 100) + (2 * 10)

        for j in xrange(tier_depth):

            data = []

            for k in xrange(len(outcomes)):

                data.append([])

                # a Counter allows us to 'reduce' all entries with similar 
                # <fp>.<alpha> pairs 
                _counter = collections.Counter()
                for key, value in _data[i][j][k]:
                    _counter[key] += value

                _counter = sorted(_counter.items())

                x = [ seq[0] for seq in _counter ]
                data[k] = [ seq[1] for seq in _counter ]

            print "X AXIS VALUES: " + str(x)

            # accumulate the values from diff. outcomes, for the same 
            # fp value 
            print data
            y_stack = np.cumsum(data, axis = 0)
            # make them represented from 0.0 to 100.0 %
            y_stack = (y_stack / max(y_stack[len(outcomes) - 1])) * 100.0
            print "Y AXIS VALUES: "
            print y_stack

            subplot_code += 1
            stackd = fig.add_subplot(subplot_code)

            stackd.fill_between(x, 0, y_stack[0,:], facecolor='lightgreen', alpha=.7)
            stackd.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor='pink', alpha=.7)
            stackd.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor='cyan', alpha=.7)
            stackd.fill_between(x, y_stack[2,:], y_stack[3,:], facecolor='gray', alpha=.7)

            stackd.set_xlabel(x_descr + str(j), fontsize=LABEL_FONT_SIZE)
            stackd.set_ylabel('Outcome occurrence %', fontsize=LABEL_FONT_SIZE)
            
            a = plt.plot([], [], color='lightgreen', linewidth=10)
            b = plt.plot([], [], color='pink', linewidth=10)
            c = plt.plot([], [], color='cyan', linewidth=10)
            d = plt.plot([], [], color='gray', linewidth=10)

            stackd.grid(True)
            fig.tight_layout()

            # x-axis handling by parts:
            stackd.set_xscale('log')
            stackd.set_xticks(x)
            stackd.set_xlim(min(x), max(x))
            # stackd.set_xticklabels(x)   
            # x_labels = stackd.get_xticklabels()
            # plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

            # y axis is easy...
            stackd.set_ylim(0.0, 100.0)

            # title for subplot
            # fig.suptitle(penalty_types[i].title())

            # legend
            stackd.legend((a[0], b[0], c[0], d[0]), ('Correct', 'Incorrect', penalty_types[i].title(), 'Dropped'), loc='lower left', fontsize=LEGEND_FONT_SIZE)

        plt.savefig(sys.argv[4] + "/stackd." + x_param_str + "." + penalty_types[i] + ".cache." + input_file_dir[-1:] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()