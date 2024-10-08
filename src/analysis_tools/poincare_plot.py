#!/usr/bin/env python
from __future__ import print_function
import sys
#from synergia.utils import Hdf5_file
import numpy
from matplotlib import pyplot
import h5py

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
    print("usage: synpoincareplot <filename1> ... <filenamen> [option1] ... [optionn] <h coord> <v coord>")
    print("available options are:")
    print("    --index=<index0>[,<index1>,...]: select indices from a bulk_tracks file (default is 0)")
    print("    --pointsize=<float>: size of plotted points (default=4.0)")
    print("    --output=<file> : save output to file (not on by default)")
    print("    --legend : give legend for each particle number in plot")

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
        if not coord in list(coords.keys()):
            do_error('Unknown coord "%s"' % coord)
    if options.indices == []:
        options.indices = [0]
    return options

def single_plot(options, particle_coords, trk):
    x = particle_coords[:, coords[options.coords[0]]]
    y = particle_coords[:, coords[options.coords[1]]]
    plot2d(x, y, options, trk)

def do_plots(options):
    # the method Figure.canvas.set_window_title() was deprecated in
    # matplotlib version 3.4 and removed in 3.6. Some github pages warn
    # that a canvas.manager might not always exist but it seems work in
    # all modern systems.
    # pyplot.figure().canvas.set_window_title('Synergia Poincare Plot')

    pyplot.figure().canvas.manager.set_window_title('Synergia Poincare Plot')
    for filename in options.inputfiles:
        #f = Hdf5_file(filename, Hdf5_file.read_only)
        f = h5py.File(filename, 'r')
        if "coords" in f.keys():
            particle_coords = f.get("coords")
            f.close()
            single_plot(options, particle_coords, 0)
        elif "track_coords" in f.keys():
            track_coords = f.get("track_coords")[()]
            #print('track_coords.shape: ', track_coords.shape)
            nturns = track_coords.shape[0]
            ntracks = track_coords.shape[1]
            mass = f.get('mass')[()]
            p_ref = f.get('pz')[()]
            #print("p_ref.shape: ", p_ref.shape)
            f.close()
            for trk in options.indices:
                particle_coords = track_coords[:,trk,0:6]
                nentries = particle_coords.shape[0]
                #print("particle_coords.shape: ", particle_coords.shape)
                #print("track_coords[:, trk, 5].shape: ", track_coords[:, trk, 5].shape)
                #print("p_ref.shape: ", p_ref.shape)
                pz = (p_ref * (1.0 + track_coords[:, trk, 5])).reshape(nentries,1)
                #print("pz.shape: ", pz.shape)
                energy = numpy.sqrt(pz*pz + mass**2).reshape(nentries,1)
                #print("energy.shape: ", energy.shape)
                #print("particle_coords.shape: ", particle_coords.shape)
                stackup = (particle_coords,pz,energy)
                #print("what is being stacked")
                #for a in stackup:
                #    print(a.shape)
                particle_coords = numpy.hstack(stackup)
                #print("particle_coords.shape: ", particle_coords.shape)
                single_plot(options, particle_coords, trk)

    pyplot.xlabel(options.coords[0])
    pyplot.ylabel(options.coords[1])
    if options.legend:
        pyplot.legend(loc='best')
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

def poincare_plot(command):
    options = handle_args(command.split())
    do_plots(options)

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_plots(options)
