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
#    bYDip=0.091*kinetic_energy/8.0
  bYDip = 0.2
  bField = numpy.array([0., bYDip, 0.], 'd')
  numTry=1000
  vv=numpy.random.rand(numTry)
  massE = 1.0e6*synergia.PH_NORM_me # Energy units here are keV..
  if (MPI.COMM_WORLD.Get_size() != 1):
    print " Not a tru MPI application!... "
    sys.exit()
  ctx=0.01 
  tInit=6.0e-9 
  for iTry in range(numTry):
     mEl.clear()
     ek= 0.01 + 1.0*vv[iTry] # in keV
     eTot = ek + massE;
     p = numpy.sqrt(eTot*eTot - massE*massE) 
     betaY = p/massE
     thisElec = numpy.array([0.0, ctx*betaY, 0.0, betaY, \
                              0.0, 1.0e-6, tInit, 0., iTry], 'd')
     if (iTry < 10): 
       print " Kinetic energy for test electron ", ek 		      
     mEl.addElectron(thisElec)
     nIter=0
     while (mEl.numLeftToTrackBC() > 0):
       mEl.propagateNoBeam(bField, "")
       if (mEl.numReachedBeamPipe() == 0):
          if (nIter == 0): 
	    print " First electron did not reached the wall... " 
          break; 
       mEl.AddFromWall(2.0, 1)
       totalKin = mEl.totalKineticEnergy()
       eLoss=ek-totalKin
       if (eLoss < -0.001): # leave room for one electron, since we at T>0... 
         print " Non energy conservation!!! Out ", totalKin, " In " , ek, " Eloss ", eLoss
	 sys.exit()
	 
       mEl.gatherNumInVaccum()
       if (iTry < 10):
         print "  After from Sec emission, between bunches iter ", \
	       nIter, " numElectrons ", mEl.numInVaccum(), \
	       " total energy ", totalKin 
       nIter+=1
       if (nIter > 25):
	 if (MPI.COMM_WORLD.Get_rank() == 0):
	   print "  Run away with one electron, cut short after 25 iteration "
	   sys.exit() 
	   break
    
     # end while on emptying buffer.. 
   
     
     
     
