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

# in honor of "stack'd", the pittsburgh restaurant next to my shadyside 
# place

LABEL_FONT_SIZE=8
LEGEND_FONT_SIZE=8

# outcome type strings
outcomes = ['Mult. hits', 'No hits', 'TP deliv', 'FP deliv.']

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    if len(sys.argv) < 2:
        print "usage: python plot-stackd-sanity-check.py <input-file-dir> <output-file-dir>"
        return

    start_time = time.time()

    input_file_dir = sys.argv[1]

    _data = []

    # we will directly fill a data array structured as 
    # _data[<outcome index>][(<scenario code>, <probability>)]
    for i in xrange(len(outcomes)):
        _data.append([])

    # read each .csv file in input_file_dir
    for file_name in os.listdir(input_file_dir):

        if file_name.endswith(".csv"):

            _file = open(input_file_dir + "/" + file_name, 'rb')

            for line in _file.readlines():
                line_splitted = line.split(",") 

                try:
                    outcome_index = int(math.log(int(line_splitted[4]), 2))
                    scenario_code = str(os.path.splitext(file_name)[0]).upper()
                    _data[outcome_index].append((scenario_code, float(line_splitted[2])))
                
                except IndexError:

                    print "line = " + line
                    print "out. index = " + str(outcome_index) + " out. str = " + outcomes[outcome_index]
                    print "scn. code = " + scenario_code
                    print "prob. = " + str(float(line_splitted[2]))

                    return

    elapsed_time = time.time() - start_time
    print "[READ FILES IN " + str(elapsed_time) + " sec]"

    fig = plt.figure()
    subplot_code = (2 * 100) + (2 * 10)

    # will save the 'stacked' data in the 'data' list, of size 
    # len(outcomes) x scenarios. the idea is for it to become somthing like:
    # 
    #         [SCN 1][SCN 2][ ... ][SCN n]
    # [OUT 1]  p11     p12    ...    p1n
    # [OUT 2]  p21     p22    ...    p2n   
    # [OUT 3]  p31     p32    ...    p3n
    # [OUT 4]  p41     p42    ...    p4n

    data = []

    for i in xrange(len(outcomes)):

        data.append([])

        # this sums all probabilities of the same <outcome, scn code> pair 
        # into a single bucket data[outcome][scn code]
        _counter = collections.Counter()
        for key, value in _data[i]:
            _counter[key] += value

        # sort [scn codes] alphabetically
        _counter = sorted(_counter.items())

        x = [ seq[0] for seq in _counter ]
        data[i] = [ seq[1] for seq in _counter ]

    print "[X AXIS VALUES]: " + str(x)

    # here, we cumulatively sum the columns of the data[] array : the idea is 
    # to sum the elements below row i to the row i, which results in an offset 
    # for the data in row i (this needs a better explanation)
    y_stack = np.cumsum(data, axis = 0)
    y_stack = (y_stack / max(y_stack[len(outcomes) - 1])) * 100.0
    print "[Y AXIS VALUES]:"
    print y_stack

    # onto the graphing...
    subplot_code += 1
    stackd = fig.add_subplot(subplot_code)

    _x = np.arange(len(x))

    stackd.fill_between(_x, 0, y_stack[0,:], facecolor='cyan', alpha=.7)
    stackd.fill_between(_x, y_stack[0,:], y_stack[1,:], facecolor='gray', alpha=.7)
    stackd.fill_between(_x, y_stack[1,:], y_stack[2,:], facecolor='lightgreen', alpha=.7)
    stackd.fill_between(_x, y_stack[2,:], y_stack[3,:], facecolor='pink', alpha=.7)

    stackd.set_xlabel("Scenarios", fontsize=LABEL_FONT_SIZE)
    stackd.set_ylabel('Outcome %', fontsize=LABEL_FONT_SIZE)
    
    a = plt.plot([], [], color='cyan', linewidth=10)
    b = plt.plot([], [], color='gray', linewidth=10)
    c = plt.plot([], [], color='lightgreen', linewidth=10)
    d = plt.plot([], [], color='pink', linewidth=10)

    stackd.grid(True)
#    fig.tight_layout()

    # x-axis handling by parts:
    stackd.set_xticks(_x)
    stackd.set_xlim(min(_x), max(_x))
    stackd.set_xticklabels(x)   
    x_labels = stackd.get_xticklabels()
    plt.setp(x_labels, rotation=90, fontsize=LABEL_FONT_SIZE)

    # y axis is easy...
    stackd.set_ylim(0.0, 100.0)

    # legend
    stackd.legend((c[0], d[0], a[0], b[0]), (outcomes[2], outcomes[3], outcomes[0], outcomes[1]), loc='lower right', fontsize=LEGEND_FONT_SIZE)

    plt.savefig(sys.argv[2] + "/sanity-check.png", bbox_inches='tight')

if __name__ == "__main__":
    main()
