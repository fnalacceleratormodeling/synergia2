#!/usr/bin/env python

import sys, os.path
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
        self.xmlfile = False
        self.line = None
        self.plots = []
        self.reader = None
        self.lfcsvfile = None
        self.lfnpfile = None
def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help(plotparams):
    print "usage: synlatticefns [options] <lattice_file> <line_name> <fn_1> ... <fn_n>"
    print "    or"
    print "       synlatticefns [options] <lattice_file.xml> <fn_1> ... <fn_n>"
    print "options governing the interpretation of the lattice are:"
    print "    --reader=[mad8|madx]"
    print "lattice from xml file is assumed by default to be a mad8 lattice.  Change"
    print "    with --reader=madx"
    print "available plots are:"
    plots = plotparams.keys()
    plots.sort()
    for plot in plots:
        print plot,
    print
    print "lattice function file save options are:"
    print "    --lfcsvfile=<filename>    (save as a csv text file)"
    print "    --lfnpfile=<filename>     (save as a numpy .npy file)"
    sys.exit(0)

def handle_args(args, plotparams):
    if len(args) < 2:
        do_help(plotparams)
    options = Options()
    args_to_remove = []
    for arg in args:
        if arg.find('--reader') == 0:
            reader = arg.split('=')[1]
            options.reader = reader
            args_to_remove.append(arg)
        elif arg.find('--lfcsvfile') ==0:
            options.lfcsvfile = arg.split("=")[1]
            args_to_remove.append(arg)
        elif arg.find('--lfnpfile') == 0:
            options.lfnpfile = arg.split("=")[1]
            args_to_remove.append(arg)
    for arg in args_to_remove:
        args.remove(arg)
    options.filename = args[0]
    # is this a lattice file or an xml file?
    if os.path.splitext(options.filename)[1] == '.xml':
        options.xmlfile = True
        options.line = "not_needed"
        start_args = 1
    else:
        options.line = args[1]
        start_args = 2
    if (len(args)<start_args+1) and (not options.lfcsvfile) and (not options.lfnpfile):
        do_help(plotparams)

    for arg in args[start_args:]:
        if arg[0] == '-':
            do_error('Unknown argument "%s"' % arg)
        else:
            if arg in plotparams.keys():
                options.plots.append(arg)
            else:
                do_error('Unknown plot "%s"' % arg)
    if options.xmlfile and not options.reader:
        options.reader = 'mad8'
        print "lattice reader defaulting to mad8"
    if not options.xmlfile and not options.reader:
        if os.path.splitext(options.filename)[1] == '.madx':
            options.reader = 'madx'
        else:
            options.reader = 'mad8'
    if (options.reader != 'madx') and (options.reader != 'mad8'):
        do_error("reader must be either 'mad8' or 'madx'")

    return options


def do_csvfile(options, lfinfo):
    lffo = open(options.lfcsvfile, "w")
    print >>lffo, "#name s alpha_x beta_x psi_x alpha_y beta_y psi_y D_x Dprime_x D_y Dprime_y"
    for lf in lfinfo:
        print >>lffo, "%16s %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g %.16g"%(lf['name'], lf['s'], lf['alpha_x'], lf['beta_x'], lf['psi_x'], lf['alpha_y'],
                                                                                    lf['beta_y'], lf['psi_y'], lf['D_x'], lf['Dprime_x'], lf['D_y'], lf['Dprime_y'])

    lffo.close()
    return


def get_lf_info():
    if options.xmlfile:
        if options.reader == 'mad8':
            adaptor_map_obj = synergia.lattice.Mad8_adaptor_map()
        elif options.reader == 'madx':
            adaptor_map_obj = synergia.lattice.MadX_adaptor_map()
        lattice = synergia.lattice.Lattice("lattice_from_xmlfile", adaptor_map_obj)
        synergia.lattice.xml_load_lattice(lattice, options.filename)
    elif options.reader == 'madx':
        lattice = synergia.lattice.MadX_reader().get_lattice(options.line, options.filename)
    else:
        lattice = synergia.lattice.Mad8_reader().get_lattice(options.line, options.filename)
    print "read lattice, ", len(lattice.get_elements()), " elements, length: ", lattice.get_length()
    lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)
    n_elem = len(lattice.get_elements())

    lf_names = ("beta_x", "alpha_x", "beta_y", "alpha_y", "psi_x", "psi_y",
                "D_x", "Dprime_x", "D_y", "Dprime_y")
    lf_info = [{}] # array of lf dicts
    index = 0

    for element in lattice.get_elements():
        index += 1
        lattice_functions = lattice_simulator.get_lattice_functions(element)
        lf = {}
        lf['name'] = element.get_name()
        lf['s'] = lattice_functions.arc_length
        lf['length'] = element.get_length()
        for lfnm in lf_names:
            lf[lfnm] = getattr(lattice_functions, lfnm)
        lf_info.append(lf)

    # set values at beginning
    lf_info[0]['s'] = 0.0
    lf_info[0]['name'] = "INITIAL"
    lf_info[0]['length'] = 0.0
    lf_info[0]['psi_x'] = 0.0
    lf_info[0]['psi_y'] = 0.0
    for lfnm in lf_names:
        lf_info[0][lfnm] = lf_info[-1][lfnm]
    return lf_info

def do_plots(options, plotparams, lf_info):
    ss = numpy.array([lfd['s'] for lfd in lf_info])
    for plot in options.plots:
        lfs = numpy.array([lfd[plotparams[plot].member] for lfd in lf_info])
        pyplot.plot(ss, lfs, label=plotparams[plot].label)
    pyplot.legend(loc='best')
    pyplot.show()

if __name__ == '__main__':
    plotparams = generate_plotparams()
    options = handle_args(sys.argv[1:], plotparams)
    lf_info = get_lf_info()
    if options.plots:
        do_plots(options, plotparams, lf_info)
    if options.lfcsvfile:
        do_csvfile(options, lf_info)
    if options.lfnpfile:
        numpy.save(options.lfnpfile, lf_info)


