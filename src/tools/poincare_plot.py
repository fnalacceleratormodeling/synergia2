#!/usr/bin/env python

import sys
import tables
from matplotlib import pyplot
import matplotlib

def plot2d(x, y, options):
    p = pyplot.plot(x, y, 'o')
    pyplot.setp(p, "markersize", float(options.point_size))
    pyplot.setp(p, "markeredgewidth", 0)


coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5

class Options:
    def __init__(self):
        self.point_size = 4.0
        self.show = True
        self.inputfiles = []
        self.outputfile = None
        self.coords = []

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synpoincareplot <filename1> ... <filenamen> [option1] ... [optionn] <h coord> <v coord>"
    print "available options are:"
    print "    --pointsize=<float>: size of plotted points (default=4.0)"
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
    if len(args) < 3:
        do_help()
    options = Options()
    first_coord = len(args) - 2
    for arg in args[:first_coord]:
        if arg[0] == '-':
            if arg == '--help':
                do_help(plotparams)
            elif arg.find('--pointsize') == 0:
                options.point_size = arg.split('=')[1]
            elif arg == '--show':
                options.show = True
            elif arg.find('--output') == 0:
                file = arg.split('=')[1]
                options.outputfile = file
                options.show = False
            else:
                do_error('Unknown option "%s"' % arg)
        else:
            options.inputfiles.append(arg)
    options.coords = args[first_coord:]
    for coord in options.coords:
        if not coord in coords.keys():
            do_error('Unknown coord "%s"' % coord)
    return options

def single_plot(options, particle_coords):
    x = particle_coords[coords[options.coords[0]], :]
    y = particle_coords[coords[options.coords[1]], :]
    plot2d(x, y, options)

def do_plots(options):
    pyplot.figure().canvas.set_window_title('Synergia Poincare Plot')
    for filename in options.inputfiles:
        f = tables.openFile(filename, 'r')
        if "coords" in dir(f.root):
            particle_coords = getattr(f.root, "coords").read()
            f.close()
            single_plot(options, particle_coords)
        elif "track_coords" in dir(f.root):
            track_coords = getattr(f.root, "track_coords").read()
            f.close()
            ntracks = track_coords.shape[0]
            for trk in range(ntracks):
                particle_coords = track_coords[trk,0:6,:]
                single_plot(options, particle_coords)

    pyplot.xlabel(options.coords[0])
    pyplot.ylabel(options.coords[1])
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_plots(options)
