#!/usr/bin/env python

import local_paths
import gourmet
import pylab
import numpy

def show_lattice_functions(lf):
    bx = pylab.plot(lf.s,lf.beta_x,'-b',linewidth=1.5)
    by = pylab.plot(lf.s,lf.beta_y,'-g',linewidth=1.5)
    pylab.legend((bx,by),('beta x','beta y'))
    pylab.xlabel('s (m)')
    pylab.ylabel('beta (m)')
    pylab.show()
    
g = gourmet.Gourmet("simplebooster.mad","cell",0.4,1.0)
g.print_elements()
lf = g.get_lattice_functions()
print "x beta function =",lf.beta_x[0],"at s =",lf.s[0]

show_lattice_functions(lf)

print "Why does this hang???"
