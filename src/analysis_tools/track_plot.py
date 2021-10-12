#!/usr/bin/env python
from __future__ import print_function
import sys
#from synergia.utils import Hdf5_file
import numpy
import h5py
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

def plot2d(x, y, coord, index):
    if index:
        label = coord + ' ' + str(index)
    else:
        label = coord
    pyplot.plot(x, y, label=label)

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5
coords['pz'] = 7
coords['energy'] = 8

class Options:
    def __init__(self):
        self.oneplot = False
        self.show = True
        self.inputfile = None
        self.outputfile = None
        self.indices = [None]
        self.coords = []

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print("usage: syntrackplot <filename> [option1] ... [optionn] <h coord1> <v coord1> ... <h coordn> <v coordn>")
    print("available options are:")
    print("    --index=<index0>[,<index1>,...]: select indices from a bulk_tracks file (default is 0)")
    print("    --oneplot : put all plots on the same axis (not on by default)")
    print("    --output=<file> : save output to file (not on by default)")
    print("    --show : show plots on screen (on by default unless --output flag is present")
    print("available coords are:")
    print("   ", end=' ')
    coord_keys = list(coords.keys())
    coord_keys.sort()
    for coord in coord_keys:
        print(coord, end=' ')
    print()
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
                do_help()
            elif arg.find('--index=') == 0:
                indices = arg.split('=')[1]
                options.indices = list(map(int,indices.split(',')))
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
            if arg in list(coords.keys()):
                options.coords.append(arg)
            else:
                do_error('Unknown coord "%s"' % arg)
    return options

#def get_particle_coords(f, options):
def get_particle_coords(h5, options):
    particle_coords = []
    if "coords" in h5.keys():
        particle_coords.append(h5.get("coords"))
    else:
        all_coords = h5.get("track_coords")
        nturns = all_coords.shape[0]
        if options.indices[0] is not None:
            indices = options.indices
        else:
            print("using default track index 0")
            indices = [0]
    mass = h5.get('mass')[()]
    p_ref = h5.get("pz")[()].reshape((nturns, 1))
    for index in indices:
        #pz = p_ref * (1.0 + all_coords[:, index, 5])
        dpop = all_coords[:, index, 5].reshape((nturns,1))
        pz = p_ref * (1.0 + dpop)
        energy = numpy.sqrt(pz*pz + mass**2)
        particle_coords.append(numpy.hstack((all_coords[:, index,:], pz, energy)))
    return particle_coords

def do_plots(options):
    #f = Hdf5_file(options.inputfile, Hdf5_file.read_only)
    h5 = h5py.File(options.inputfile, 'r')
    rows, cols = get_layout(len(options.coords))
    pyplot.figure().canvas.set_window_title('Synergia Track Viewer')
    #all_particle_coords = get_particle_coords(f, options)
    all_particle_coords = get_particle_coords(h5, options)
    for particle_coords, index in zip(all_particle_coords, options.indices):
        plot_index = 1
        for coord in options.coords:
            #x = f.read_array1d("s")
            x = h5.get("s")[()]
            y = particle_coords[:,coords[coord]][()]
            if not options.oneplot:
                pyplot.subplot(rows, cols, plot_index)
            plot2d(x, y, coord, index)
            plot_index += 1
            pyplot.legend()
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_plots(options)
