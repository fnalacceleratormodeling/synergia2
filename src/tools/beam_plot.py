#!/usr/bin/env python

import sys
import tables
import numpy
import math
import pylab
from matplotlib import pyplot
from mpl_toolkits.axes_grid import make_axes_locatable

def plot_density(x, y, label, bins):
    fancylabel = label.replace('_', ' ')

    fig = pyplot.figure(1)
    axScatter = pyplot.subplot(111)
    divider = make_axes_locatable(axScatter)
    axHistx = divider.new_vertical(1.2, pad=0.1, sharex=axScatter)
    axHisty = divider.new_horizontal(1.2, pad=0.1, sharey=axScatter)
    fig.add_axes(axHistx)
    fig.add_axes(axHisty)

    pyplot.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
         visible=False)

    H, xedges, yedges = numpy.histogram2d(x, y, bins=bins)
    for i in range(0, H.shape[0]):
        for j in range(0, H.shape[1]):
            if H[i, j] > 0:
                H[i, j] = math.log(H[i, j])
    X, Y = pylab.meshgrid(xedges, yedges)
    axScatter.pcolor(X, Y, H.transpose())
    axScatter.set_axis_bgcolor((0, 0, 0.5))
    axHistx.hist(x, bins=bins)
    axHisty.hist(y, bins=bins, orientation='horizontal')
#    axScatter.plot(x, y, '.', label=fancylabel)

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
        self.show = True
        self.hcoord = None
        self.vcoord = None
        self.bins = 10

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synbeamplot <filename> [option1] ... [optionn] <h coord> <v coord>"
    print "available options are:"
    print "    --nohist : do not show histograms (not on by default)"
    print "    --bins=<num> : number of bins in each direction"
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
            elif arg.find('--bins') == 0:
                options.bins = int(arg.split('=')[1])
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
    plot_density(particles[:, coords[options.hcoord]],
                 particles[:, coords[options.vcoord]], 'foobar', options.bins)
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

if __name__ == '__main__':
#    plotparams = generate_plotparams()
    options = handle_args(sys.argv[1:])
    do_plots(options)
