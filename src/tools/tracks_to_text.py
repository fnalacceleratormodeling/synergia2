#!/usr/bin/env python

import sys, os.path
import tables

class Options:
    def __init__(self):
        self.header = 2
        self.filename = None

def do_error(message):
    sys.stderr.write(message + '\n')
    sys.exit(1)

def do_help():
    print "usage: syntrack2txt [option] <filename>"
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
    turns = ifile.root.repetition.read()
    npart = ifile.root.track_coords.shape[0]
    nturns = ifile.root.repetition.shape[0]

    ofilename = os.path.splitext(options.filename)[0] + '.txt'
    ofile = open(ofilename, 'w')
    if options.header > 1:
        ofile.write('''# units:
# x, y, cdt : meters
# xp = px/pref
# yp = py/pref
# dp = (delta ptotal)/pref
''')
    if options.header > 0:
        ofile.write('#')
        ofile.write('%10s'%'turn')
	for p in range(npart):
            for col in ['x', 'xp', 'y', 'py', 'cdt', 'dp', 'id']:
                ofile.write('%24s' % col)
        ofile.write('\n')

    if ifile.root.track_coords.shape[2] != nturns:
        raise RuntimeError, "number of turns doesn't agree between arrays repetition and track_coords"

    for t in range(len(turns)):
        ofile.write(' ')
        ofile.write('%10d'%turns[t])
        for part in range(npart):
            for coord in range(7):
                ofile.write('%24.16e' % ifile.root.track_coords[part, coord, t])
        ofile.write('\n')
    ofile.close()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_conversion(options)
