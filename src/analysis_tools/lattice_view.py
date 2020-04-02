#!/usr/bin/env python

import sys, os.path
import synergia
from matplotlib import pyplot
from math import sin, cos, pi
import numpy

xmlfile = False
used_legends = []
def plot_element(element, x, y, angle, attributes, highlight):
    global used_legends
    l = element.get_length()
    bend_angle = element.get_bend_angle()
    if element.get_length() == 0.0:
        xnew = x
        ynew = y
        anglenew = angle
    else:
        ancestors = element.get_ancestors()
        type = element.get_type()
        if type in attributes:
            attribute = attributes[type]
        else:
            attribute = attributes['default']
        color = list(attribute.color)
        if (ancestors.count(highlight) == 0) and (highlight != None):
            offset = 0.75
            color[0] = offset + (1.0 - offset) * color[0]
            color[1] = offset + (1.0 - offset) * color[1]
            color[2] = offset + (1.0 - offset) * color[2]
        label = None
        if used_legends.count(type) == 0:
            label = type
            used_legends.append(type)
        if element.get_bend_angle() == 0.0:
            xnew = x + l * cos(angle)
            ynew = y + l * sin(angle)
            anglenew = angle
            if label:
                pyplot.plot([x, xnew], [y, ynew], '-',
                            color=color,
                            label=label,
                            linewidth=attribute.width)
            else:
                pyplot.plot([x, xnew], [y, ynew], '-',
                            color=color,
                            linewidth=attribute.width)
        else:
            num = 8
            xn = numpy.zeros([num + 1], numpy.float64)
            yn = numpy.zeros([num + 1], numpy.float64)
            xn[0] = x
            yn[0] = y
            for i in range(1, num + 1):
                angle += bend_angle / num
                xn[i] = xn[i - 1] + (l / num) * cos(angle)
                yn[i] = yn[i - 1] + (l / num) * sin(angle)
            xnew = xn[num]
            ynew = yn[num]
            anglenew = angle
            if label:
                pyplot.plot(xn, yn, '-',
                            color=color,
                            label=label,
                            linewidth=attribute.width)
            else:
                pyplot.plot(xn, yn, '-',
                            color=color,
                            linewidth=attribute.width)

    return xnew, ynew, anglenew

class Attributes:
    def __init__(self, color, width=2.0):
        self.color = color
        self.width = width

def get_attributes():
    attributes = {}
    attributes['drift'] = Attributes((0.5, 0.5, 0.5))
    attributes['quadrupole'] = Attributes((1, 0, 0))
    attributes['sbend'] = Attributes((0, 0, 1))
    attributes['rbend'] = Attributes((0, 0, 0.5))
    attributes['sextupole'] = Attributes((0, 0.5, 0.5))
    attributes['octupole'] = Attributes((1, 0.5, 0))
    attributes['monitor'] = Attributes((1, 1, 0))
    attributes['hmonitor'] = Attributes((0.5, 1, 0))
    attributes['vmonitor'] = Attributes((1, 0.5, 0))
    attributes['kicker'] = Attributes((0, 1, 1))
    attributes['hkicker'] = Attributes((0, 1, 0.5))
    attributes['vkicker'] = Attributes((0, 0.5, 1))
    attributes['rfcavity'] = Attributes((0.5,0,0.5))
    attributes['default'] = Attributes((0,0,0),1.0)
    return attributes

class Options:
    def __init__(self):
        self.lattice_file = None
        self.lattice = None
        self.legend = True
        self.highlight = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print("usage: synlatticeview <filename> <lattice> [option1] ... [optionn]")
    print("available options are:")
    print("    --nolegend : do not show legend")
    print("    --highlight=<name> : highlight line <name>")
    print("    --reader=<type> : use lattice reader <type> (\"mad8\" or \"madx\")")
    print("                    : default reader is madx for *.madx files, mad8 otherwise")
    sys.exit(0)

def handle_args(args):
    if len(args) < 1:
        do_help()
    options = Options()
    options.lattice_file = args[0]
    # is this an xml file?
    if os.path.splitext(options.lattice_file)[1] == '.xml':
        options.xmlfile = True
        start_args = 1
        # default reader mad8 for xml files can be changed with --reader=
        options.reader = 'mad8'
    else:
        if len(args) < 2:
            do_help()
        options.xmlfile = False
        start_args = 2
        options.lattice = args[1]
    options.reader = None
    for arg in args[start_args:]:
        if arg == '--help':
            do_help(plotparams)
        elif arg == '--nolegend':
            options.legend = False
        elif arg.find('--highlight') == 0:
            highlight = arg.split('=')[1]
            options.highlight = highlight
        elif arg.find('--reader') == 0:
            reader = arg.split('=')[1]
            options.reader = reader
        else:
            do_error('Unknown argument "%s"' % arg)
    if not options.reader:
        if os.path.splitext(options.lattice_file)[1] == '.madx':
            options.reader = 'madx'
        else:
            options.reader = 'mad8'
    if (options.reader != 'madx') and (options.reader != 'mad8'):
        do_error("reader must be either 'mad8' or 'madx'")
    return options

def do_plot(options):
    if options.xmlfile:
        lattice = synergia.lattice.Lattice()
        synergia.lattice.xml_load_lattice(lattice, options.lattice_file)
    elif not options.xmlfile and options.reader == 'madx':
        lattice = synergia.lattice.MadX_reader().get_lattice(options.lattice,
                                                             options.lattice_file)
    else:
        lattice = synergia.lattice.Mad8_reader().get_lattice(options.lattice,
                                                             options.lattice_file)
    attributes = get_attributes()
    pyplot.figure().canvas.set_window_title('Synergia Lattice Viewer')
    x = 0.0
    y = 0.0
    angle = 0.0
    count = 0
    for element in lattice.get_elements():
        x, y, angle = plot_element(element, x, y, angle, attributes,
                                   options.highlight)
        count += 1
    print("total number of elements =", count)
    print("total angle =", lattice.get_total_angle(),"rad (%g degrees)" % \
        (lattice.get_total_angle()*180.0/pi))

    pyplot.axes().set_aspect('equal', 'datalim')
    try:
        pyplot.axes().margins(0.1, 0.1)
    except:
        # margins not available in older versions of Matplotlib
        pass
    if options.legend:
        pyplot.legend()
    pyplot.show()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_plot(options)
