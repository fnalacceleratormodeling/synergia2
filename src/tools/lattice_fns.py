#!/usr/bin/env python

import sys, os.path
import tables
from matplotlib import pyplot
import synergia
import numpy

def plot2d(x, y, label):
    fancylabel = label.replace('_', ' ')
    pyplot.plot(x, y, label=fancylabel)

class Params:
    def __init__(self, member, label):
        self.member = member
        self.label = label

def generate_plotparams():
    plotparams = {}
    plotparams['betax'] = Params('beta_x', r'$\beta_x$')
    plotparams['betay'] = Params('beta_y', r'$\beta_y$')
    plotparams['alphax'] = Params('alpha_x', r'$\alpha_x$')
    plotparams['alphay'] = Params('alpha_y', r'$\alpha_y$')
    plotparams['psix'] = Params('psi_x', r'$\psi_x$')
    plotparams['psiy'] = Params('psi_y', r'$\psi_y$')
    plotparams['Dx'] = Params('D_x', r'$D_x$')
    plotparams['Dy'] = Params('D_y', r'$D_y$')
    plotparams['Dprimex'] = Params('Dprime_x', r'$D^{\prime}_x$')
    plotparams['Dprimey'] = Params('Dprime_y', r'$D^{\prime}_y$')
    return plotparams

class Options:
    def __init__(self):
        self.filename = None
        self.line = None
        self.plots = []
        self.reader = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help(plotparams):
    print "usage: synlatticefns <lattice_file> <line_name> <fn_1> ... <fn_n>"
    print "available plots are:"
    print "   ",
    plots = plotparams.keys()
    plots.sort()
    for plot in plots:
        print plot,
    print
    sys.exit(0)

def handle_args(args, plotparams):
    if len(args) < 3:
        do_help(plotparams)
    options = Options()
    for arg in args:
        if arg.find('--reader') == 0:
            reader = arg.split('=')[1]
            options.reader = reader
            args.remove(arg)
    options.filename = args[0]
    options.line = args[1]
    for arg in args[2:]:
        if arg[0] == '-':
            do_error('Unknown argument "%s"' % arg)
        else:
            if arg in plotparams.keys():
                options.plots.append(arg)
            else:
                do_error('Unknown plot "%s"' % arg)
    if not options.reader:
        if os.path.splitext(options.filename)[1] == '.madx':
            options.reader = 'madx'
        else:
            options.reader = 'mad8'
    if (options.reader != 'madx') and (options.reader != 'mad8'):
        do_error("reader must be either 'mad8' or 'madx'")

    return options

def do_plots(options, plotparams):
    if options.reader == 'madx':
        lattice = synergia.lattice.MadX_reader().get_lattice(options.line, options.filename)
    else:
        lattice = synergia.lattice.Mad8_reader().get_lattice(options.line, options.filename)
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)
    for plot in options.plots:
        n = len(lattice.get_elements())
        ss = numpy.zeros([n+1], numpy.float64)
        lfs = numpy.zeros([n+1], numpy.float64)
        lattice_functions = None
        index = 0
        for element in lattice.get_elements():
            index += 1
            lattice_functions = lattice_simulator.get_lattice_functions(element)
            ss[index] = lattice_functions.arc_length
            lfs[index] = getattr(lattice_functions, plotparams[plot].member)
        ss[0] = 0.0
        if (plot == 'psix') or (plot == 'psiy'):
            lfs[0] = 0.0
        else:
            lfs[0] = getattr(lattice_functions, plotparams[plot].member)
        pyplot.plot(ss, lfs, label=plotparams[plot].label)
    pyplot.legend(loc='best')
    pyplot.show()

if __name__ == '__main__':
    plotparams = generate_plotparams()
    options = handle_args(sys.argv[1:], plotparams)
    do_plots(options, plotparams)
