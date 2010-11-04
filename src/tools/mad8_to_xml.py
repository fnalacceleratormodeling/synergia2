#!/usr/bin/env python
import sys
import synergia

if len(sys.argv) != 4:
    print "usage:",
    print sys.argv[0], "<mad8 file> <line name> <xml file>"
    print "    Reads line <line name> from <mad8 file> and writes to <xml file>."
    sys.exit(1)

mad8_file = sys.argv[1]
line = sys.argv[2]
xml_file = sys.argv[3]

reader = synergia.lattice.Mad8_reader()
reader.parse(mad8_file)
print "found lines:", reader.get_lines()
lattice = reader.get_lattice(line, enable_cache_read=False,
                             enable_cache_write=False)
synergia.lattice.xml_save_lattice(lattice, xml_file)
print "wrote", line, "to", xml_file
