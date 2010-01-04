#!/usr/bin/env python

import local_paths
import gourmet
import pylab
import numpy

g = gourmet.Gourmet("simplebooster.mad","cell",0.4)
g.insert_space_charge_markers(4)
g.generate_maps(200.0e6)
i = 0
for map in g.maps:
    i += 1
    print "map %d:"%i
    print numpy.array2string(map,precision=4,suppress_small=1)
    print

print "Why does this hang???"
