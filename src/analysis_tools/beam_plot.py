#!/usr/bin/env python

from __future__ import print_function
import sys
import numpy
import h5py
import math
import pylab
from matplotlib import pyplot
from mpl_toolkits.axes_grid1.axes_divider import make_axes_locatable

def plot_density(x, y, label, options):
    fancylabel = label.replace('_', ' ')

    fig,axScatter  = pyplot.subplots(1,1)
    fig.suptitle('Synergia Phase Space Distribution')
    divider = make_axes_locatable(axScatter)
    axHistx = divider.new_vertical(1.2, pad=0.1, sharex=axScatter)
    axHisty = divider.new_horizontal(1.2, pad=0.1, sharey=axScatter)
    fig.add_axes(axHistx)
    fig.add_axes(axHisty)

    pyplot.setp(axHistx.get_xticklabels() + axHisty.get_yticklabels(),
         visible=False)

    H, xedges, yedges = numpy.histogram2d(x, y, bins=options.bins)
    for i in range(0, H.shape[0]):
        for j in range(0, H.shape[1]):
            if H[i, j] > 0:
                H[i, j] = math.log(H[i, j])
    # contour plots use the center of the bins, not the edges
    if options.contour:
        xedges = (xedges[0:-1] + xedges[1:])/2.0
        yedges = (yedges[0:-1] + yedges[1:])/2.0
    X, Y = pylab.meshgrid(xedges, yedges)
    if not options.contour:
        axScatter.pcolor(X, Y, H.transpose())
    else:
        if options.num_contour:
            # select number of contours
            axScatter.contourf(X, Y, H.transpose(), options.num_contour)
        else:
            # matplotlib chooses number of contours
            axScatter.contourf(X, Y, H.transpose())

    axScatter.set_facecolor((0, 0, 0.5)) if hasattr(axScatter, 'set_facecolor') else axScatter.set_axis_bgcolor((0, 0, 0.5))
    axHistx.hist(x, bins=options.bins)
    axHisty.hist(y, bins=options.bins, orientation='horizontal')
#    axScatter.plot(x, y, '.', label=fancylabel)

    # if limits were specified, set the axis limits to match
    xlims = list(axScatter.get_xlim())
    ylims = list(axScatter.get_ylim())
    if options.minh != -sys.float_info.max:
        xlims[0] = options.minh
    if options.maxh != sys.float_info.max:
        xlims[1] = options.maxh
    if options.minv != -sys.float_info.max:
        ylims[0] = options.minv
    if options.maxv != sys.float_info.max:
        ylims[1] = options.maxv
    axScatter.set_xlim(xlims)
    axScatter.set_ylim(ylims)
    axScatter.set_xlabel(options.hcoord, fontdict={'fontsize': 'large'})
    axScatter.set_ylabel(options.vcoord, fontdict={'fontsize': 'large'})

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['cdt'] = 4
coords['dpop'] = 5
coords['pz'] = 7
coords['energy'] = 8
coords['t'] = 9
coords['z'] = 10

class Options:
    def __init__(self):
        self.hist = True
        self.inputfile = None
        self.outputfile = None
        self.show = True
        self.hcoord = None
        self.vcoord = None
        self.bins = 10
        self.minh = -sys.float_info.max
        self.maxh = sys.float_info.max
        self.minv = -sys.float_info.max
        self.maxv = sys.float_info.max
        self.contour = None
        self.num_contour = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print("usage: synbeamplot <filename> [option1] ... [optionn] <h coord> <v coord>")
    print("available options are:")
    print("    --nohist : do not show histograms (not on by default)")
    print("    --bins=<num> : number of bins in each direction")
    print("    --minh=<num> : minimum limit on horizontal axis data")
    print("    --maxh=<num> : maximum limit on horizontal axis data")
    print("    --minv=<num> : minimum limit on vertical axis data")
    print("    --maxv=<num> : maximum limit on vertical axis data")
    print("    --contour    : draw contours instead of color density")
    print("    --contour=<num>: draw <num> levels of contours")
    print("    --output=<file> : save output to file (not on by default)")
    print("    --show : show plots on screen (on by default unless --output flag is present")
    print("available coords are:")
    print("   ", end=' ')
    coord_list = list(coords.keys())
    coord_list.sort()
    for coord in coord_list:
        print(coord, end=' ')
    print()
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
                options.hist = False
            elif arg.find('--bins') == 0:
                options.bins = int(arg.split('=')[1])
            elif arg == '--show':
                options.show = True
            elif arg.find('--output') == 0:
                file = arg.split('=')[1]
                options.outputfile = file
                options.show = False
            elif arg.find('--minh') == 0:
                options.minh = float(arg.split('=')[1])
            elif arg.find('--maxh') == 0:
                options.maxh = float(arg.split('=')[1])
            elif arg.find('--minv') == 0:
                options.minv = float(arg.split('=')[1])
            elif arg.find('--maxv') == 0:
                options.maxv = float(arg.split('=')[1])
            elif arg.find('--contour') == 0:
                options.contour = True
                if arg.find('=') >= 0:
                    options.num_contour = int(arg.split('=')[1])
            else:
                do_error('Unknown argument "%s"' % arg)
    return options

def do_plots(options):
    c = 299792458.0
    h5 = h5py.File(options.inputfile, 'r')
    #f = Hdf5_file(options.inputfile, 'r')
    particles = h5.get('particles')
    masks = h5.get('particles_masks')
    npart = particles.shape[0]
    mass = h5.get('mass')[()]
    p_ref = h5.get('pz')[()]
    betagamma = p_ref/mass
    gamma = numpy.sqrt(betagamma**2 + 1)
    beta = betagamma/gamma
    #print "p_ref: ", p_ref
    #print "mass: ", mass
    pz = (p_ref * (1.0 + particles[:,5])).reshape(npart, 1)
    #print "pz.shape: ", pz.shape
    energy = numpy.sqrt(pz*pz + mass**2).reshape(npart, 1)
    time = (particles[:,4]*1.0e9/c).reshape(npart,1)
    #print "energy.shape: ", energy.shape
    z = (particles[:,4]*beta).reshape(npart,1)
    particles = numpy.hstack((particles, pz, energy,time, z))
    
    selected_particles = ( (masks[:] != 0) * (particles[:, coords[options.hcoord]] >= options.minh) *
                       (particles[:, coords[options.hcoord]] < options.maxh) *
                       (particles[:, coords[options.vcoord]] >= options.minv) *
                       (particles[:, coords[options.vcoord]] < options.maxv))

    plot_density(particles[selected_particles, coords[options.hcoord]],
                 particles[selected_particles, coords[options.vcoord]], 'foobar', options)
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

if __name__ == '__main__':
#    plotparams = generate_plotparams()
    options = handle_args(sys.argv[1:])
    do_plots(options)
