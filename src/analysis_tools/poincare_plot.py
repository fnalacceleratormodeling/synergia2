#!/usr/bin/env python

import sys
from synergia.utils import Hdf5_file
import numpy
from matplotlib import pyplot

def plot2d(x, y, options, trk):
    if options.legend:
        p = pyplot.plot(x, y, 'o', label="%d"%trk)
    else:
        p = pyplot.plot(x, y, 'o')
    pyplot.setp(p, "markersize", float(options.point_size))
    pyplot.setp(p, "markeredgewidth", 0.1)


coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5
coords['pz'] = 6
coords['energy'] = 7

class Options:
    def __init__(self):
        self.point_size = 4.0
        self.show = True
        self.inputfiles = []
        self.outputfile = None
        self.coords = []
        self.indices = []
        self.legend = False

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synpoincareplot <filename1> ... <filenamen> [option1] ... [optionn] <h coord> <v coord>"
    print "available options are:"
    print "    --index=<index0>[,<index1>,...]: select indices from a bulk_tracks file (default is 0)"
    print "    --pointsize=<float>: size of plotted points (default=4.0)"
    print "    --output=<file> : save output to file (not on by default)"
    print "    --legend : give legend for each particle number in plot"

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
                do_help()
            elif arg.find('--index') == 0:
                indices = arg.split('=')[1]
                for index in indices.split(','):
                    options.indices.append(int(index))
            elif arg.find('--pointsize') == 0:
                options.point_size = arg.split('=')[1]
            elif arg == '--legend':
                options.legend = True
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
    if options.indices == []:
        options.indices = [0]
    return options

def single_plot(options, particle_coords, trk):
    x = particle_coords[:, coords[options.coords[0]]]
    y = particle_coords[:, coords[options.coords[1]]]
    plot2d(x, y, options, trk)

def do_plots(options):
    pyplot.figure().canvas.set_window_title('Synergia Poincare Plot')
    for filename in options.inputfiles:
        f = Hdf5_file(filename, Hdf5_file.read_only)
        if "coords" in f.get_member_names():
            particle_coords = f.read_array2d("coords")
            f.close()
            single_plot(options, particle_coords, 0)
        elif "track_coords" in f.get_member_names():
            track_coords = f.read_array3d("track_coords")
            nturns = track_coords.shape[0]
            ntracks = track_coords.shape[1]
            mass = f.read_double('mass')
            p_ref = f.read_array1d('pz').reshape(nturns,1)
            f.close()
            for trk in options.indices:
                particle_coords = track_coords[:,trk,0:6]
                pz = p_ref * (1.0 + track_coords[:, trk, 5]).reshape(nturns,1)
                #print "pz.shape: ", pz.shape
                energy = numpy.sqrt(pz*pz + mass**2)
                #print "energy.shape: ", energy.shape
                particle_coords = numpy.hstack((particle_coords,pz,energy))
                #print "particle_coords.shape: ", particle_coords.shape
                single_plot(options, particle_coords, trk)

    pyplot.xlabel(options.coords[0])
    pyplot.ylabel(options.coords[1])
    if options.legend:
        pyplot.legend(loc='best')
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_plots(options)
