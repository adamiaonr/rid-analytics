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

def main():

    start_time = time.time()

    if len(sys.argv) < 2:
        print "usage: python plot-outcome-vs-dist.py <input-scn-file> <input-file-dir> <output-file-dir>"
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

            tier_depth = int(splitted[1])

        elif splitted[0] == 'fp_prob':

            for v in splitted[1].split(","):

                fp_prob.append(float(v.rstrip()))

        elif splitted[0] == 'alpha':

            for v in splitted[1].split(","):

                alpha.append(float(v.rstrip()))

    # prepare the arrays which will be used as placeholders for data
    # we will extract all data into a multidimensional array:
    #   - [penalty][fp prob][alpha][outcome] -> (dist, probability)
    _data = []

    for i in xrange(len(penalty_types)):
        _data.append([])

        for j in xrange(len(fp_prob)):
            _data[i].append([])

            for k in xrange(len(alpha)):
                _data[i][j].append([])

                for l in xrange(len(outcomes)):
                    _data[i][j][k].append([])

    # dir where input data .csv files are
    input_file_dir = sys.argv[2]

    # read the outcomes .csv files for each cache distance
    for i in xrange(tier_depth):        

        # build outcomes .csv file name
        outcomes_file_name = input_file_dir + "/cache." + str(i + 1) + "/outcomes.csv"

        if not os.path.isfile(outcomes_file_name):
            print "ERROR: " + outcomes_file_name + " doesn't exist. abort."
            return

        # open the file for a particular outcome & penalty combination
        outcomes_file = open(outcomes_file_name, 'rb')

        for line in outcomes_file.readlines():
            
            splitted = line.split(",")

            try:
                fp_prob.index(float(splitted[0]))
                alpha.index(float(splitted[1]))

            except ValueError:
                continue

            try:
                # remember: the order is:
                #   - [penalty][fp prob][alpha][outcome] -> (dist, probability)

                _data[penalty_types.index(str(splitted[2]))][fp_prob.index(float(splitted[0]))][alpha.index(float(splitted[1]))][outcomes.index(str(splitted[4]))].append((int(i + 1), float(splitted[5])))

            # catch any index out of bounds error        
            except IndexError:

                try:

                    print "{0:penalty} :"   + str(penalty_types.index(str(splitted[2])))
                    print "{1:fp prob} :"   + str(fp_prob.index(float(splitted[0])))
                    print "{2:alpha} :"     + str(alpha.index(float(splitted[1])))
                    print "{3:outcome} :"   + str(outcomes.index(str(splitted[4])))

                except ValueError:

                    print "{0:penalty} :"   + splitted[2]
                    print "{1:fp prob} :"   + splitted[0]
                    print "{2:alpha} :"     + splitted[1]
                    print "{3:outcome} :"   + splitted[4]

                return

    elapsed_time = time.time() - start_time

    print "[READ FILE(S) IN " + str(elapsed_time) + " sec]"

    # done with the file reading. draw stack'd graphs for each penalty case.
    for i in xrange(len(penalty_types)):

        fig = plt.figure()
        subplot_code = (2 * 100) + (2 * 10)

        for j in xrange(len(fp_prob)):

            for k in xrange(len(alpha)):

                data = []

                for l in xrange(len(outcomes)):

                    data.append([])

                    # a Counter allows us to 'reduce' all entries with similar 
                    # <fp>.<alpha> pairs 
                    _counter = collections.Counter()
                    for key, value in _data[i][j][k][l]:
                        _counter[key] += value

                    _counter = sorted(_counter.items())

                    x = [ seq[0] for seq in _counter ]
                    data[l] = [ seq[1] for seq in _counter ]

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
                chart = fig.add_subplot(subplot_code)

                if k > 0:

                    chart.fill_between(x, 0, y_stack[0,:], facecolor='green', alpha=.5)
                    chart.fill_between(x, y_stack[0,:], y_stack[1,:], facecolor='red', alpha=.5)
                    chart.fill_between(x, y_stack[1,:], y_stack[2,:], facecolor='blue', alpha=.5)
                    chart.fill_between(x, y_stack[2,:], y_stack[3,:], facecolor='gray', alpha=.5)

                    chart.set_xlabel("Cache distance (in hops)", fontsize=LABEL_FONT_SIZE)
                    chart.set_ylabel('Outcome occurrence %', fontsize=LABEL_FONT_SIZE)

                    # y axis is easy...
                    chart.set_ylim(0.0, 100.0)

                    a = plt.plot([], [], color='green', linewidth=5, alpha=0.5)
                    b = plt.plot([], [], color='red', linewidth=5, alpha=0.5)
                    c = plt.plot([], [], color='blue', linewidth=5, alpha=0.5)

                    chart.legend((a[0], b[0], c[0]), ('Correct', 'Incorrect', penalty_types[i].title()), loc='upper left', fontsize=LEGEND_FONT_SIZE)

                else:

                    chart.plot(x, data[0] / np.cumsum(data, axis = 0)[len(outcomes) - 1], color='green', linewidth=2.5, alpha=0.5)
                    chart.plot(x, data[1] / np.cumsum(data, axis = 0)[len(outcomes) - 1], color='red', linewidth=2.5, alpha=0.5)
                    chart.plot(x, data[2] / np.cumsum(data, axis = 0)[len(outcomes) - 1], color='blue', linewidth=2.5, alpha=0.5)

                    chart.set_xlabel("Cache distance (in hops)", fontsize=LABEL_FONT_SIZE)
                    chart.set_ylabel('Prob. of occurrence [0.0, 1.0]', fontsize=LABEL_FONT_SIZE)

                    chart.set_yscale('log')
                    #print math.log10(min(min(d[1:]) for d in np.cumsum(data, axis = 0)))
                    #chart.set_ylim(math.pow(10, math.floor(math.log10(min(min(d[1:]) for d in np.cumsum(data, axis = 0))))), 1.0)
                    chart.set_ylim(math.pow(10, -16), 1.0)
                
                    a = plt.plot([], [], color='green', linewidth=5, alpha=0.5)
                    b = plt.plot([], [], color='red', linewidth=5, alpha=0.5)
                    c = plt.plot([], [], color='blue', linewidth=5, alpha=0.5)
                    #d = plt.plot([], [], color='gray', linewidth=10)

                    chart.legend((a[0], b[0], c[0]), ('Correct', 'Incorrect', penalty_types[i].title()), loc='upper right', fontsize=LEGEND_FONT_SIZE)

                chart.grid(True)
                fig.tight_layout()

                # x-axis handling by parts:
                chart.set_xticks(x)
                chart.set_xlim(min(x), max(x))

                # for _x in x:
                #     x_labels_str.append('%-.1e' % _x)
                chart.set_xticklabels(x, fontsize=LABEL_FONT_SIZE)   
                #x_labels = chart.get_xticklabels()
                #plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

                chart.set_title("FP rate " + str(fp_prob[j]) + " / Alpha " + str(alpha[k]), fontsize=LABEL_FONT_SIZE)

                # legend
                #chart.legend((a[0], b[0], c[0], d[0]), ('Correct', 'Incorrect', penalty_types[i].title(), 'Dropped'), loc='lower left', fontsize=LEGEND_FONT_SIZE)

        plt.savefig(sys.argv[3] + "/outcome-vs-dist." + penalty_types[i] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()