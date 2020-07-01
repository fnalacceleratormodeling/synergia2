#!/usr/bin/env python

import sys, os
import synergia

if len(sys.argv) != 3:
    print "usage: ",sys.argv[0]," <lattice_file>  <line/sequence name>"

lattice_file = sys.argv[1]
line_name = sys.argv[2]

(lat_pre, lat_ext) = os.path.splitext(lattice_file)

if (lat_ext == ".lat") or (lat_ext == ".mad"):
    print "reading line ", line_name, " from ", lattice_file
    # this is a mad8 lattice file
    lattice = synergia.lattice.Mad8_reader().get_lattice(line_name, lattice_file)
    print "Read lattice ", len(lattice.get_elements()), " elements"
elif (lat_ext == ".madx") or (lat_ext == ".seq"):
    print "reading line/sequence ", line_name, " from ", lattice_file
    # this is a madX lattice or sequence file
    lattice = synergia.lattice.MadX_reader().get_lattice(line_name, lattice_file)
    print "Read lattice ", len(lattice.get_elements()), " elements"

lsxfile = lat_pre+".lsx"
synergia.utils.write_lsexpr_file(lattice.as_lsexpr(), lsxfile)
print "Write lattice to file ", lsxfile
