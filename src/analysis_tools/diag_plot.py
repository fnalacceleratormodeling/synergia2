#!/usr/bin/env python

import sys
from matplotlib import pyplot
from synergia.utils import Hdf5_file

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

def plot2d(x, y, label, extra_label):
    fancylabel = label.replace('_', ' ')
    if extra_label:
        fancylabel += extra_label
    pyplot.plot(x, y, label=fancylabel)

class Params:
    def __init__(self, label, y_attr, y_index1=None, y_index2=None):
        self.label = label
        self.y_attr = y_attr
        self.y_index1 = y_index1
        self.y_index2 = y_index2

coords = {}
coords['x'] = 0
coords['xp'] = 1
coords['y'] = 2
coords['yp'] = 3
coords['z'] = 4
coords['zp'] = 5

def generate_plotparams():
    plotparams = {}
    for label in coords.keys():
        for label2 in coords.keys():
            if coords[label2] > coords[label]:
                corr = label + '_' + label2 + '_corr'
                plotparams[corr] = Params(corr, 'corr',
                                           coords[label], coords[label2])
                mom2 = label + '_' + label2 + '_mom2'
                plotparams[mom2] = Params(mom2, 'mom2',
                                           coords[label], coords[label2])
        std = label + '_std'
        plotparams[std] = Params(std, 'std', coords[label])
        mean = label + '_mean'
        plotparams[mean] = Params(mean, 'mean', coords[label])
    plotparams['x_emit'] = Params('x_emit', 'emitx')
    plotparams['y_emit'] = Params('y_emit', 'emity')
    plotparams['z_emit'] = Params('z_emit', 'emitz')
    plotparams['xy_emit'] = Params('xy_emit', 'emitxy')
    plotparams['xyz_emit'] = Params('xyz_emit', 'emitxyz')
    plotparams['particles'] = Params('particles', 'num_particles')
    return plotparams

class Options:
    def __init__(self):
        self.ind_var = 's'
        self.oneplot = False
        self.show = True
        self.inputfiles = []
        self.outputfile = None
        self.legends = True
        self.plots = []

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help(plotparams):
    print "usage: syndiagplot <filename> [option1] ... [optionn] <plot1> ... <plotn>"
    print "available options are:"
    print "    --userep: use repetition instead of s for independent variable"
    print "    --oneplot : put all plots on the same axis (not on by default)"
    print "    --nolegend : suppress legends (not on by default)"
    print "    --output=<file> : save output to file (not on by default)"
    print "    --show : show plots on screen (on by default unless --output flag is present"
    print "available plots are:"
    print "   ",
    plots = plotparams.keys()
    plots.sort()
    for plot in plots:
        print plot,
    print
    sys.exit(0)

def handle_args(args, plotparams):
    if len(args) < 2:
        do_help(plotparams)
    options = Options()
    filenames = args[0]
    options.inputfiles = filenames.split(',')
    for arg in args[1:]:
        if arg[0] == '-':
            if arg == '--help':
                do_help(plotparams)
            elif arg == '--userep':
                options.ind_var = 'repetition'
            elif arg == '--oneplot':
                options.oneplot = True
            elif arg == '--nolegends':
                options.legends = False
            elif arg == '--show':
                options.show = True
            elif arg.find('--output') == 0:
                file = arg.split('=')[1]
                options.outputfile = file
                options.show = False
            else:
                do_error('Unknown argument "%s"' % arg)
        else:
            if arg in plotparams.keys():
                options.plots.append(arg)
            else:
                do_error('Unknown plot "%s"' % arg)
    return options

def hdf5_read_any(hdf5_file, member):
    try:
        the_type = hdf5_file.get_atomic_type(member)
    except:
        sys.stderr.write('syndiagplot: data member "%s" not found in Hdf5 file\n'
                        % member)
        sys.exit(1)

    dims = hdf5_file.get_dims(member)
    if the_type == Hdf5_file.double_type:
        if len(dims) == 0:
            return hdf5_file.read_double(member)
        elif len(dims) == 1:
            return hdf5_file.read_array1d(member)
        elif len(dims) == 2:
            return hdf5_file.read_array2d(member)
        elif len(dims) == 3:
            return hdf5_file.read_array3d(member)
    elif the_type == Hdf5_file.int_type:
        if len(dims) == 0:
            return hdf5_file.read_int(member)
        elif len(dims) == 1:
            return hdf5_file.read_array1i(member)

def do_plot(inputfile, options, plotparams, multiple_files):
    f = Hdf5_file(inputfile, Hdf5_file.read_only)
    rows, cols = get_layout(len(options.plots))
    plot_index = 1
    y_label = ""
    for plot in options.plots:
        params = plotparams[plot]
        x = hdf5_read_any(f, options.ind_var)
        ymaster = hdf5_read_any(f, params.y_attr)
        if (params.y_index1 == None) and (params.y_index2 == None):
            y = ymaster
        elif (params.y_index2 == None):
            y = ymaster[params.y_index1, :]
        else:
            y = ymaster[params.y_index1, params.y_index2, :]
        if not options.oneplot:
            pyplot.subplot(rows, cols, plot_index)
        extra_label = None
        if multiple_files:
            extra_label = ' ' + inputfile
        plot2d(x, y, plot, extra_label)
        plot_index += 1
        pyplot.xlabel(options.ind_var)
        if options.oneplot:
            if y_label == "":
                y_label = plot
            else:
                y_label += "," + plot
        else:
            y_label = plot
        pyplot.ylabel(y_label)
        if options.legends:
            pyplot.legend(loc='best')

def do_plots(options, plotparams):
    pyplot.figure().canvas.set_window_title('Synergia Diagnostics')
    multiple_files = False
    if len(options.inputfiles) > 1:
        multiple_files = True
    for inputfile in options.inputfiles:
        do_plot(inputfile, options, plotparams, multiple_files)
    if options.outputfile:
        pyplot.savefig(options.outputfile)
    if options.show:
        pyplot.show()

if __name__ == '__main__':
    plotparams = generate_plotparams()
    options = handle_args(sys.argv[1:], plotparams)
    do_plots(options, plotparams)
