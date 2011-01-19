#!/usr/bin/env python

import sys
import tables
from matplotlib import pyplot

def get_layout(num):
    if num == 1:
        return 1, 1
    elif num == 2:
        return 2, 1
    elif num == 3:
        return 3, 1
    elif num == 4:
        return 2, 2
    elif num <= 6:
        return 2, 2
    elif num <= 9:
        return 3, 3
    elif num <= 12:
        return 3, 4
    elif num <= 16:
        return 4, 4
    else:
        do_error("Too many plots")

def plot2d(x, y, label):
    pyplot.plot(x, y, 'o', label=label)

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5

class Options:
    def __init__(self):
        self.oneplot = False
        self.show = True
        self.inputfile = None
        self.outputfile = None
        self.coords = []

    def get_pairs(self):
        pairs = []
        for i in range(0, len(self.coords), 2):
            pair = [self.coords[i], self.coords[i + 1]]
            pairs.append(pair)
        return pairs

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synpoincareplot <filename> [option1] ... [optionn] <h coord1> <v coord1> ... <h coordn> <v coordn>"
    print "available options are:"
    print "    --oneplot : put all plots on the same axis (not on by default)"
    print "    --output=<file> : save output to file (not on by default)"
    print "    --show : show plots on screen (on by default unless --output flag is present"
    print "available coords are:"
    print "   ",
    coord_keys = coords.keys()
    coord_keys.sort()
    for coord in coord_keys:
        print coord,
    print
    sys.exit(0)

def handle_args(args):
    if len(args) < 2:
        do_help()
    options = Options()
    filename = args[0]
    options.inputfile = filename
    for arg in args[1:]:
        if arg[0] == '-':
            if arg == '--help':
                do_help(plotparams)
            elif arg == '--oneplot':
                options.oneplot = True
            elif arg == '--show':
                options.show = True
            elif arg.find('--output') == 0:
                file = arg.split('=')[1]
                options.outputfile = file
                options.show = False
            else:
                do_error('Unknown argument "%s"' % arg)
        else:
            if arg in coords.keys():
                options.coords.append(arg)
            else:
                do_error('Unknown coord "%s"' % arg)
    if len(options.coords) % 2 == 1:
        do_error("coordinates must be specified in pairs")
    return options

def do_plots(options):
    f = tables.openFile(options.inputfile, 'r')
    rows, cols = get_layout(len(options.get_pairs()))
    plot_index = 1
    particle_coords = getattr(f.root, "coords").read()
    for pair in options.get_pairs():
        x = particle_coords[coords[pair[0]], :]
        y = particle_coords[coords[pair[1]], :]
        if not options.oneplot:
            pyplot.subplot(rows, cols, plot_index)
        plot2d(x, y, pair[0] + "," + pair[1])
        plot_index += 1
        pyplot.legend()
    f.close()
    pyplot.show()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_plots(options)
