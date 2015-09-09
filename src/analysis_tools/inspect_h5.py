#!/usr/bin/env python

import sys
from synergia.utils import Hdf5_file

def usage(error=False):
    print "usage:"
    print "     syninspecth5 <hdf5_file>"
    print "     " + "    to list members"
    print "or"
    print "     syninspecth5 <hdf5_file> <member>"
    print "     " + "      to display member"
    if error:
        retval = 1
    else:
        retval = 0
    sys.exit(retval)

def print_entries(filename):
    f = Hdf5_file(filename, Hdf5_file.read_only)
    for entry in  f.get_member_names():
        if entry[0] != '_':
            print entry

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

def print_entry(filename, entry):
    f = Hdf5_file(filename, Hdf5_file.read_only)
    value = hdf5_read_any(f, entry)
    print entry,
    if hasattr(value,"shape"):
        print "shape =", value.shape,
    print
    print value

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
