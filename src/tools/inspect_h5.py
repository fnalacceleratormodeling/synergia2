#!/usr/bin/env python

import sys
import tables

def usage(error=False):
    print "usage:"
    print "     syninspecth5 <hdf5_file>"
    print "     " + "    to list members"
    print "or"
    print "     syninspecth5 <hdf5_file> <entry>"
    print "     " + "      to display entry"
    if error:
        retval = 1
    else:
        retval = 0
    sys.exit(retval)

def print_entries(filename):
    f = tables.openFile(filename, 'r')
    for entry in  dir(f.root):
        if entry[0] != '_':
            print entry
    f.close()

def print_entry(filename, entry):
    f = tables.openFile(filename, 'r')
    value =  getattr(f.root, entry).read()
    print entry,
    if hasattr(value,"shape"):
        print "shape =", value.shape,
    print
    print value
    f.close()

if __name__ == '__main__':
    if len(sys.argv) < 2:
        usage(error=True)
    elif len(sys.argv) == 2:
        if sys.argv[1] == '--help':
            usage()
        else:
            print_entries(sys.argv[1])
    elif len(sys.argv) == 3:
        print_entry(sys.argv[1], sys.argv[2])
    else:
        usage(error=True)
