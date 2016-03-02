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
        print "usage: python plot-outcome-vs-alpha.py <input-scn-file> <input-file-dir> <output-file-dir>"
        return

    # read .scn file
    scn_file = open(sys.argv[1], 'rb')

    # extract the tier_depth, alpha and fp_prob values from the .scn file
    tier_depth = 0

    fp_prob = []
    alpha = []

    for line in scn_file.readlines():

        # extract the var names (before the '=' sign)
        splitted = line.split("=")

        # extract the var values
        if splitted[0] == 'tier_depth':

            tier_depth = splitted[1]

        elif splitted[0] == 'fp_prob':

            for v in splitted[1].split(","):

                fp_prob.append(float(v.rstrip()))

        elif splitted[0] == 'alpha':

            for v in splitted[1].split(","):

                alpha.append(float(v.rstrip()))

    # dir where input data .csv files are
    input_file_dir = sys.argv[2]

    # build outcomes .csv file name
    outcomes_file_name = input_file_dir + "/outcomes.csv"

    if not os.path.isfile(outcomes_file_name):
        print "ERROR: " + outcomes_file_name + " doesn't exist. abort."
        return

    # prepare the arrays which will be used as placeholders for data
    # we will extract all data into a multidimensional array:
    #   - [penalty][fp prob][outcome] -> (alpha, probability)
    _data = []

    for i in xrange(len(penalty_types)):
        _data.append([])

        for j in xrange(len(fp_prob)):
            _data[i].append([])

            for k in xrange(len(outcomes)):
                _data[i][j].append([])

    # open the file for a particular outcome & penalty combination
    outcomes_file = open(outcomes_file_name, 'rb')

    prob_sum = 0.0

    for line in outcomes_file.readlines():
        
        splitted = line.split(",")

        try:
            # remember: the order is:
            #   - [penalty][fp prob][outcome] -> (alpha, probability)
            _data[penalty_types.index(str(splitted[2]))][fp_prob.index(float(splitted[0]))][outcomes.index(str(splitted[4]))].append((float(splitted[1]), float(splitted[5])))

        # catch any index out of bounds error        
        except IndexError:

            try:

                print "{0:penalty} :"   + str(penalty_types.index(str(splitted[2])))
                print "{1:fp prob} :"   + str(fp_prob.index(float(splitted[0])))
                print "{2:outcome} :"   + str(outcomes.index(str(splitted[4])))

            except ValueError:

                print "{0:penalty} :"   + splitted[2]
                print "{1:fp prob} :"   + splitted[0]
                print "{2:outcome} :"   + splitted[4]

            return

    elapsed_time = time.time() - start_time

    print "[READ FILE IN " + str(elapsed_time) + " sec]"

    # done with the file reading. draw stack'd graphs for each penalty case.
    for i in xrange(len(penalty_types)):

        fig = plt.figure()
        subplot_code = (2 * 100) + (2 * 10)

        for j in xrange(len(fp_prob)):

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

            stackd.fill_between(x, 0, y_stack[0,:], facecolor='green', alpha=.5)
            stackd.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor='red', alpha=.5)
            stackd.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor='blue', alpha=.5)
            stackd.fill_between(x, y_stack[2,:], y_stack[3,:], facecolor='gray', alpha=.5)

            # stackd.fill_between(x, 0, y_stack[0,:], facecolor='lightgreen', alpha=.7)
            # stackd.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor='pink', alpha=.7)
            # stackd.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor='cyan', alpha=.7)
            # stackd.fill_between(x, y_stack[2,:], y_stack[3,:], facecolor='red', alpha=.7)

            stackd.set_xlabel("Alpha", fontsize=LABEL_FONT_SIZE)
            stackd.set_ylabel('Outcome occurrence %', fontsize=LABEL_FONT_SIZE)
            
            a = plt.plot([], [], color='green', linewidth=10, alpha=.5)
            b = plt.plot([], [], color='red',   linewidth=10, alpha=.5)
            c = plt.plot([], [], color='blue',  linewidth=10, alpha=.5)
#            d = plt.plot([], [], color='gray', linewidth=10)

            stackd.grid(True)
            fig.tight_layout()

            # x-axis handling by parts:
            stackd.set_xscale('log')
            stackd.set_xticks(x)
            stackd.set_xlim(min(x), max(x))

            # for _x in x:
            #     x_labels_str.append('%-.1e' % _x)

            stackd.set_xticklabels(x, fontsize=LABEL_FONT_SIZE)   
            #x_labels = stackd.get_xticklabels()
            #plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

            # y axis is easy...
            stackd.set_ylim(0.0, 100.0)

            stackd.set_title("FP Rate " + str(fp_prob[j]), fontsize=LABEL_FONT_SIZE)

            # legend
#            stackd.legend((a[0], b[0], c[0], d[0]), ('Correct', 'Incorrect', penalty_types[i].title(), 'Dropped'), loc='lower left', fontsize=LEGEND_FONT_SIZE)
            stackd.legend((a[0], b[0], c[0]), ('Correct', 'Incorrect', penalty_types[i].title()), loc='upper left', fontsize=LEGEND_FONT_SIZE)

        plt.savefig(sys.argv[3] + "/outcome-vs-alpha." + penalty_types[i] + ".cache." + input_file_dir[-1:] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()