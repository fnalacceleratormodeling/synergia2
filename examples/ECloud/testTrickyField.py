import numpy
import time
import math
import os
import sys

import synergia
from mpi4py import MPI
import txphysics.txionpack
from electronFlock import ElectronFlock


# Testing of the generation on the wall... 

if ( __name__ == '__main__'):
  
  mEl = ElectronFlock()
  mEl.setMaxXDimPipe(0.055)
  mEl.setMaxYDimPipe(0.025)
  mEl.setNumberTrajectoryDump(10)
#    bYDip=0.091*kinetic_energy/8.0
  
  bYDip = 0.2
  bField = numpy.array([0., bYDip, 0.], 'd')
  mEl.setMagnetModel(4, 20.)
  numTry=6
  massE = 1.0e6*synergia.PH_NORM_me # Energy units here are keV..
  if (MPI.size != 1):
    print " Not a tru MPI application!... "
    sys.exit()
  ctx=0.01 
  tInit=6.0e-9
  
  for iTry in range(numTry):
     mEl.clear()
     tokenTraj="testTrickyFieldb_" + "%d_" % iTry
     ek= 0.015
     eTot = ek + massE;
     p = numpy.sqrt(eTot*eTot - massE*massE) 
     betaY = p/massE
     zz=-0.3 + iTry*0.1
     thisElec = numpy.array([0.01, ctx*betaY, 0.01, betaY, \
                              zz, 0.3*betaY, tInit, 0., iTry], 'd')
     mEl.addElectron(thisElec)
     mEl.propagateNoBeam(bField, tokenTraj)
     # end while on emptying buffer.. 
   
     
     
     
