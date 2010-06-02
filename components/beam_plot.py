#!/usr/bin/env python

import sys
import tables
from matplotlib import pyplot

def plot_density(x, y, label):
    fancylabel = label.replace('_', ' ')
    pyplot.plot(x, y, '.', label=fancylabel)

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5

class Options:
    def __init__(self):
        self.hist = True
        self.inputfile = None
        self.outputfile = None
        self.hcoord = None
        self.vcoord = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: beam_plot.py <filename> [option1] ... [optionn] <h coord> <v coord>"
    print "available options are:"
    print "    --nohist : do not show histograms (not on by default)"
    print "    --output=<file> : save output to file (not on by default)"
    print "    --show : show plots on screen (on by default unless --output flag is present"
    print "available coords are:"
    print "   ",
    coord_list = coords.keys()
    coord_list.sort()
    for coord in coord_list:
        print coord,
    print
    sys.exit(0)

def handle_args(args):
    if len(args) < 2:
        do_help()
    options = Options()
    filename = args[0]
    options.inputfile = filename
    options.hcoord = args[len(args) - 2]
    options.vcoord = args[len(args) - 1]
    for arg in args[1:len(args) - 2]:
        if arg[0] == '-':
            if arg == '--help':
                do_help()
            elif arg == '--nohist':
                options.hist = false
            elif arg == '--show':
                options.show = True
            elif arg.find('--output') == 0:
                file = arg.split('=')[1]
                options.outputfile = file
                options.show = False
            else:
                do_error('Unknown argument "%s"' % arg)
    return options

def do_plots(options):
    f = tables.openFile(options.inputfile, 'r')
    particles = f.root.particles.read()
    plot_density(particles[:, coords[options.hcoord]], particles[:, coords[options.vcoord]],'foobar')
    pyplot.show()

if __name__ == '__main__':
#    plotparams = generate_plotparams()
    options = handle_args(sys.argv[1:])
    do_plots(options)
