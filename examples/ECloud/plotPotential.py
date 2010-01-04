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
  phi.read_from_file("Potential_Grid_64_Nodes_8.dat")
  numGC=phi.get_points().get_length() # I think!.... 
  totalQ=2.0e11*proton_charge/numGC  
  plotPotential.plotPotentialX(phi, totalQ, "G64N8")
#
  phi=Real_scalar_field()
  phi.read_from_file("Potential_Grid_64_Nodes_1.dat")
  numGC=phi.get_points().get_length() # I think!.... 
  totalQ=2.0e11*proton_charge/numGC  
  plotPotential.plotPotentialX(phi, totalQ, "G64N1")
 
def plotPotentialX(phi, factor, token):
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
  dx=1.05*sizes[0]/n
  dy=1.05*sizes[1]/n
  dz=0.8*sizes[2]/n
  localFact=(8.98755e9/0.04954)*factor
  x = -0.45*sizes[0]
  y = -0.45*sizes[1]
#  z = -0.87*sizes[2]
  z=0.
  while (i < n):
    rr[i]=numpy.sqrt(x*x + y*y + z*z)
    loc = numpy.array([x, y, z],'d')
    vals[i]= localFact*phi.get_val(loc)
    ders[i]= localFact*phi.get_deriv(loc, 0)
    i=i+1
    x+=dx
    y+=dy
#    z+=dz
  pylab.plot(1000.0*rr, vals, 'b-', label='Phi')
  pylab.plot(1000.0*rr, -1.0*ders/100., 'r--',label='-Ex/100' )
  pylab.xlabel('x (mm)')
  pylab.ylabel('phi (Volt) or E (V/m)')
  pylab.legend()
  fName="Potential"+token+".pdf"
  pylab.savefig(fName)
  pylab.clf()
#  pylab.show()
def plotPotentialVertE(phi, factor, token):
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
  dx=1.05*sizes[0]/n
  dy=1.05*sizes[1]/n
  dz=0.8*sizes[2]/n
  localFact=(8.98755e9/0.04954)*factor
  y = -0.45*sizes[1]
#  z = -0.87*sizes[2]
  z=0.
  x= 0.015
  while (i < n):
    loc = numpy.array([x, y, z],'d')
    rr[i] = y;
    vals[i]= localFact*phi.get_val(loc)
    ders[i]= localFact*phi.get_deriv(loc, 0)
    i=i+1
#    x+=dx
    y+=dy
#    z+=dz
  pylab.plot(1000.0*rr, vals, 'b-', label='Phi')
  pylab.plot(1000.0*rr, -1.0*ders/100., 'r--',label='-Ex/100' )
  pylab.xlabel('x (mm)')
  pylab.ylabel('phi (Volt) or E (V/m)')
  pylab.legend()
  fName="VertEFieldProf"+token+".pdf"
  pylab.savefig(fName)
  pylab.clf()
#  pylab.show()
   
      
