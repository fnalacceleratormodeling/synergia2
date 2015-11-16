#!/usr/bin/env python
import sys
import synergia

if len(sys.argv) != 4:
    print "usage:",
    print "synmad8tolsx <mad8 file> <line name> <lsx file>"
    print "    Reads line <line name> from <mad8 file> and writes to <lsx file>."
    sys.exit(1)

mad8_file = sys.argv[1]
line = sys.argv[2]
lsx_file = sys.argv[3]

reader = synergia.lattice.Mad8_reader()
reader.parse(mad8_file)
print "found lines:", reader.get_lines()
lattice = reader.get_lattice(line, enable_cache_read=False,
                             enable_cache_write=False)
synergia.utils.write_lsexpr_file(lattice.as_lsexpr(), lsx_file)
print "wrote", line, "to", lsx_file
