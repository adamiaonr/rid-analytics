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

def custom_ceil(x, base=5):
#    return int(base * math.ceil(float(x)/base))
    return math.ceil(float(x)) + 1

def main():

    log_scale = 1
    markers = ['b-', 'r-', 'c-', 'm-', 'k-', 'y-', 'g-']

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

    r_size = 7 
    c_size = 6

    fp_rates = [[0 for x in range(c_size)] for x in range(r_size)] 

    i = 0
    j = 0
    labels = []

    for line in f.readlines():

        j = 0
        splitted = line.split(",")

        for split in splitted: 

            if j == 0:
                labels.append(split)
            else:
                fp_rates[i][j - 1] = float(split)

            j += 1

        i += 1

    print fp_rates
    print labels

    # let's plot the damn thing then...
    fig = plt.figure()
    fprate_plot = fig.add_subplot(111)
    fprate_plot.grid(True)

    y_max = 1.0
    y_min = y_max

    lines = [0 for x in range(7)]
    print lines

    print range(c_size)

    for i in xrange(7):

        _y_min = math.pow(10, math.floor(math.log10(min(fp_rates[i]))))
        if _y_min < y_min:
            y_min = _y_min

        fprate_plot.plot([1, 2], [fp_rates[i][0], fp_rates[i][0]], markers[i], linewidth=1.5, label=labels[i])

        ar = np.array(range(c_size))
        fprate_plot.plot(ar + 2, fp_rates[i], markers[i], linewidth=1.5)
        fprate_plot.plot([c_size + 1, c_size + 2], [fp_rates[i][c_size - 1], fp_rates[i][c_size - 1]], markers[i], linewidth=1.5)
        fprate_plot.plot(ar + c_size + 2, fp_rates[i][::-1], markers[i], linewidth=1.5)
        fprate_plot.plot([2 * c_size + 1, 2 * c_size + 2], [fp_rates[i][0], fp_rates[i][0]], markers[i], linewidth=1.5)

    if (log_scale):
        fprate_plot.set_yscale('log')
    else:
        y_min = 0.0

    fprate_plot.set_xlabel('Hop', fontsize=LABEL_FONT_SIZE)
    fprate_plot.set_ylabel('FP rate', fontsize=LABEL_FONT_SIZE)
    fprate_plot.legend(loc='lower right', fontsize=LEGEND_FONT_SIZE)
    
    fprate_plot.set_xlim(1.0, 2 * c_size + 2)
    fprate_plot.set_xticks(np.arange(1.0, 2 * c_size + 3))
    fprate_plot.set_ylim(y_min, y_max)

#    plt.show()
    plt.savefig(sys.argv[2], bbox_inches='tight')


if __name__ == "__main__":
    main()
 

