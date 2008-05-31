#
# Plot the potential
#
import numpy
import time
import math
import os
import sys
import pylab

import synergia
import s2_fish
from s2_containers import *
import plotPotential
#
proton_charge = 1.60217646e-19

if ( __name__ == '__main__'):
  phi=Real_scalar_field()
  phi.read_from_file("Potential2e11_A1.dat")
  numGC=64*64*128 # I think!.... 
  totalQ=2.0e11*proton_charge/numGC  
  plotPotential.plotPotentialX(phi, totalQ)

 
def plotPotentialX(phi, factor):
          # Check the arguments.  Expect a Real scalar field 
  argName1=str(phi.__class__)
  if (argName1.find("Real_scalar_field") == -1): 
    print " plotPotential: Wrong 1st argument type, expect a Real_scalar_field, type is  ", \
	    argName1
    print " Fatal Error "     
    sys.exit() 

  i=0
  n=500
  r=0.
  vals=numpy.zeros(n)
  ders=numpy.zeros(n)
  rr=numpy.zeros(n)
  sizes=phi.get_physical_size()
  dr=0.4*sizes[0]/n
  localFact=(8.98755e9/0.04954)*factor
  while (i < n):
    rr[i]=r
    loc = numpy.array([rr[i], rr[i], 0.],'d')
    vals[i]= localFact*phi.get_val(loc)
    ders[i]= localFact*phi.get_deriv(loc, 0)
    i=i+1
    r=r+dr
  pylab.plot(1000.0*rr, vals, 'b-', label='Phi')
  pylab.plot(1000.0*rr, -1.0*ders/100., 'r--',label='-Ex/100' )
  pylab.xlabel('x (mm)')
  pylab.ylabel('phi (Volt) or E (V/m)')
  pylab.legend()
  pylab.savefig('PotentialProp0v2.pdf')
  pylab.clf()
#  pylab.show()
   
      
