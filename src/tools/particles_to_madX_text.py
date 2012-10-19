#!/usr/bin/env python

import sys, os.path
import tables
import numpy as np

class Options:
    def __init__(self):
        self.header = 2
        self.filename = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: synpart2txt [option] <filename>"
    print "available options are:"
    print "    --long-header : add long header describing units and column labels"
    print "    --short-header : add short header of column labels"
    print "    --no-header : just include coordinates"
    sys.exit(0)

def handle_args(args):
    if len(args) < 1:
        do_help()
    options = Options()
    filename = args[0]
    options.inputfile = filename
    options.hcoord = args[len(args) - 2]
    options.vcoord = args[len(args) - 1]
    for arg in args:
        if arg[0] == '-':
            if arg == '--help':
                do_help()
            elif arg == '--long-header':
                options.header = 2
            elif arg == '--short-header':
                options.header = 1
            elif arg == '--no-header':
                options.header = 0
            else:
                do_error('Unknown argument "%s"' % arg)
        else:
            if options.filename == None:
                options.filename = arg
            else:
                do_error('Extra filename argument "%s"' % arg)
    if options.filename == None:
        do_error("No filename specified")
    return options

def do_conversion(options):
    ifile = tables.openFile(options.filename, 'r')
    ofilename = os.path.splitext(options.filename)[0] + '-madX.txt'
    ofile = open(ofilename, 'w')
    if options.header > 1:
        ofile.write('''# units:
# madX canonical variables from
# http://mad.web.cern.ch/mad/madx.old/Introduction/tables.html#canon
# X, Y, T : meters
# PX = px/pref
# PY = py/pref
# T = -c (delta-T)
# PT = (delta E)/pref
''')
    m = ifile.root.mass[()]
    pz = ifile.root.pz[()]
    E0 = np.sqrt(pz**2 + m**2)
    m_over_p0 = m/pz
    E0_over_p0 = E0/pz

    particles = ifile.root.particles.read()

    ifile.close()

    # convert synergia to MAD-X
    # c*dt -> -c*dt
    # delta-P/P -> delta-E/p
    particles[:,4] = -particles[:,4]
    # sqrt( (dP/P + 1)**2 + (m/P)**2 ) - E0/p
    tarray = particles[:,5] + 1.0
    particles[:,5] = np.sqrt(tarray*tarray + m_over_p0**2) - E0_over_p0

    ofile.write('# mass = %g [GeV/c^2], pref = %g [GeV/c]\n' % (m,pz))
    if options.header > 0:
        ofile.write('#')
        for col in ['X', 'PX', 'Y', 'PY', ' T ', 'PT', 'id']:
            ofile.write('%24s' % col)
        ofile.write('\n')

    for pnum in range(particles.shape[0]):
        part = particles[pnum,:]
        ofile.write(' ')
        for coord in part:
            ofile.write('%24.16e' % coord)
        ofile.write('\n')
    ofile.close()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_conversion(options)
