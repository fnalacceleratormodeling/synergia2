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
    ofilename = os.path.splitext(options.filename)[0] + '.txt'
    ofile = open(ofilename, 'w')
    if options.header > 1:
        ofile.write('''# units:
# x, y, cdt : meters
# xp = px/pref
# yp = py/pref
# dp = (delta ptotal)/pref
''')
        pz = ifile.root.pz.read()
        ofile.write('# pref = %g [GeV/c]\n' % pz)
    if options.header > 0:
        ofile.write('#')
        for col in ['x', 'xp', 'y', 'py', 'cdt', 'dp', 'id']:
            ofile.write('%24s' % col)
        ofile.write('\n')
    for part in ifile.root.particles.iterrows():
        ofile.write(' ')
        for coord in part:
            ofile.write('%24.16e' % coord)
        ofile.write('\n')
    ofile.close()

if __name__ == '__main__':
    options = handle_args(sys.argv[1:])
    do_conversion(options)
