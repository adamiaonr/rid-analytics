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
LEGEND_FONT_SIZE=8

outcomes = ['correct dest.', 'wrong dest. > orig. server', 'fp detect. > orig. server', 'dropped']
penalty_types = ['feedback', 'fallback']

def main():

    start_time = time.time()

    if len(sys.argv) < 2:
        print "usage: python plot-outcome-vs-hop.py <input-scn-file> <input-file-dir> <output-file-dir>"
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
    #   - [penalty][fp prob][alpha][outcome] -> (max. tier, probability)
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

    # build outcomes .csv file name
    outcomes_file_name = input_file_dir + "/outcomes.csv"

    if not os.path.isfile(outcomes_file_name):
        print "ERROR: " + outcomes_file_name + " doesn't exist. abort."
        return

    # open the file for a particular outcome & penalty combination
    outcomes_file = open(outcomes_file_name, 'rb')

    max_hops = 0

    for line in outcomes_file.readlines():
        
        splitted = line.split(",")

        try:
            fp_prob.index(float(splitted[0]))
            alpha.index(float(splitted[1]))

        except ValueError:
            continue

        try:
            # remember: the order is:
            #   - [penalty][fp prob][alpha][outcome] -> (max. tier, probability)

            _data[penalty_types.index(str(splitted[2]))][fp_prob.index(float(splitted[0]))][alpha.index(float(splitted[1]))][outcomes.index(str(splitted[4]))].append((int(splitted[6]), float(splitted[5])))

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
                prob_sum = 0.0

                for l in xrange(len(outcomes)):

                    data.append([])

                    # a Counter allows us to 'reduce' all entries with similar 
                    # <fp>.<alpha> pairs 
                    _counter = collections.Counter()
                    for key, value in _data[i][j][k][l]:
                        prob_sum += value
                        _counter[key] += value

                    print "[PROB SUM = " + str(prob_sum) + " ]"

                    _counter = sorted(_counter.items())

                    x = [ seq[0] for seq in _counter ]
                    data[l] = [ seq[1] for seq in _counter ]

                    print "X AXIS VALUES: " + str(x)

                # outcome probabilities for the same number of hops
                print data

                subplot_code += 1
                chart = fig.add_subplot(subplot_code)

                _x = np.arange(1.0, len(x) + 1)

                chart.bar(_x - 0.375,  (data[0] / np.sum(data[0], axis = 0)), 0.25, color='green',   alpha=0.5)
                chart.bar(_x - 0.125,  (data[1] / np.sum(data[1], axis = 0)), 0.25, color='red',     alpha=0.5)
                chart.bar(_x + 0.125,  (data[2] / np.sum(data[2], axis = 0)), 0.25, color='blue',    alpha=0.5)

                # chart.plot(x, (data[0] / np.sum(data[0])) * 100.0, color='green',   linewidth=2.5, alpha=0.5)
                # chart.plot(x, (data[1] / np.sum(data[1])) * 100.0, color='red',     linewidth=2.5, alpha=0.5)
                # chart.plot(x, (data[2] / np.sum(data[2])) * 100.0, color='blue',    linewidth=2.5, alpha=0.5)

                chart.set_xlabel("Max. tier", fontsize=LABEL_FONT_SIZE)
                chart.set_ylabel("Probability [0, 1]", fontsize=LABEL_FONT_SIZE)

                #print math.log10(min(min(d[1:]) for d in np.cumsum(data, axis = 0)))
                #chart.set_ylim(math.pow(10, math.floor(math.log10(min(min(d[1:]) for d in np.cumsum(data, axis = 0))))), 1.0)
                #chart.set_ylim(0.0, math.ceil(max(max(d[1:]) for d in data) * 10.0) / 10.0 + 0.1)
                chart.set_ylim(0.0, 1.0)
            
                a = plt.plot([], [], color='green', linewidth=5, alpha=0.5)
                b = plt.plot([], [], color='red', linewidth=5, alpha=0.5)
                c = plt.plot([], [], color='blue', linewidth=5, alpha=0.5)
                #d = plt.plot([], [], color='gray', linewidth=10)

                chart.legend((a[0], b[0], c[0]), ('Correct', 'Incorrect', penalty_types[i].title()), loc='upper center', fontsize=LEGEND_FONT_SIZE)

                chart.grid(True)
                fig.tight_layout()

                # x-axis handling by parts:
                chart.set_xticks(_x)
                chart.set_xlim(0.0, max(_x) + 1)

                # for _x in x:
                #     x_labels_str.append('%-.1e' % _x)
                chart.set_xticklabels(x, fontsize=LABEL_FONT_SIZE)   
                #x_labels = chart.get_xticklabels()
                #plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

                chart.set_title("FP rate " + str(fp_prob[j]) + " / Alpha " + str(alpha[k]), fontsize=LABEL_FONT_SIZE)

                # legend
                #chart.legend((a[0], b[0], c[0], d[0]), ('Correct', 'Incorrect', penalty_types[i].title(), 'Dropped'), loc='lower left', fontsize=LEGEND_FONT_SIZE)

        plt.savefig(sys.argv[3] + "/outcome-vs-tier." + penalty_types[i] + ".cache." + input_file_dir[-1:] + ".png", bbox_inches='tight')

if __name__ == "__main__":
    main()