#!/usr/bin/env python
# Electron flock, 
import numpy
import synergia
import s2_fish
import ECloudPy
import sys

class ElectronFlock:
#
    def __init__(self):
    
      self.data = [] # intented to be a list of numpy arrays, 7 doubles.
                     # 6 for phase space (3*(x,beta_x), time)
      self.dataBP = [] # Electrons stuck on the beam pipe,7 doubles.
      self.nBad=0  # Number of electron lost in inf. loop or other propagation errors. 
      self.mass = 1.0e6*synergia.PH_NORM_me # Energy units here are keV..
      self.useRK = 1
      self.timeSliceWb=0.2 # don't recall what that is about..
      self.totalCharge=0. # in Coulomb
#      
#   Initialize the GSL random number for Laundau distribution.
#
      ECloudPy.initGslRan()

    def setTotalCharge(self, charge):
      self.totalCharge=charge
      
    def addElectron(self, anElectron):
    # No check of argument, my excuse, no time for this.. 
      self.data.append(anElectron)

    def clear(self):
      self.data = []
   
    # Adding electron.  Over one meter.. 
    def addFromGas(self, numIons, prescaleFact, betaProtons, \
                   xOffset,  yOffset, xWidth, yWidth, bunchLength): 
      gamProtons = numpy.sqrt(1.0/(1.0-betaProtons*betaProtons))
      massElecOMassProton = synergia.PH_NORM_me/synergia.PH_NORM_mp
      	   
# prescale factor of 100, for now...     
      num_electrons =  int(numIons*prescaleFact)
    
      aGaussX=numpy.random.normal(xOffset, xWidth, 100)
      aGaussY=numpy.random.normal(yOffset, yWidth, 100)
      aGaussZ=numpy.random.normal(0., bunchLength, 100)
      EmaxElec = 2.0* self.mass * 1.0e3 * (betaProtons * betaProtons/(gamProtons*gamProtons))/ \
                (1.0 + 2.0* gamProtons * massElecOMassProton + \
		massElecOMassProton*massElecOMassProton) # in eV
      for i in range(num_electrons):
        ii = i%100  # by chunks of 100 ... Could go faster ? 
        if ((ii == 0) & (i > 1)):
      	  aGaussX=numpy.random.normal(xOffset, xWidth, 100)
          aGaussY=numpy.random.normal(yOffset, yWidth, 100)
          aGaussZ=numpy.random.normal(0., bunchLength, 100)
	 # Assume for now average energy of ~ 3 eV, Landau distributed.
        
        eElectKin = ECloudPy.getLandauEnergyDist(3., EmaxElec) # in eV
        eElect = eElectKin + 1.0e9*synergia.PH_NORM_me # in eV
        betaElec = numpy.sqrt(2.*eElectKin*1.0e-9/synergia.PH_NORM_me)	 
        if (eElectKin > 100.):  # relativistic case.. 
	   gamElecInv = 1.0e9*synergia.PH_NORM_me/eElect
	   betaElec = numpy.sqrt(1.0 - gamElecInv*gamElecInv)
	 # isotropic    
        dct = 2.0*numpy.random.rand() - 1.
        betaz = betaElec*dct
        sct = numpy.sqrt(1.0-dct*dct)
        dphi = 2.0*numpy.pi*numpy.random.rand()
        betax = betaElec*sct*numpy.sin(dphi) 
        betay = betaElec*sct*numpy.cos(dphi) 
        thisElec = numpy.array([aGaussX[ii], betax, aGaussY[ii], betay, \
                              aGaussZ[ii], betaz, 0.])
#      print "an Elec", thisElec
        self.addElectron(thisElec)
    
      self.setTotalCharge(synergia.PH_MKS_e * numIons)
    
    def numInVaccum(self):
      return len(self.data) 
      
    def numReachedBeamPipe(self):
      return len(self.dataBP) 
      
    def averageKineticEnergy(self):  # in keV.. 
      bb=0.
      if (len(self.data) < 2): 
        return 0.
      for ee in self.data:
        bb+=numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3] + ee[5]*ee[5])
      bb/=len(self.data)
      gam = 1./numpy.sqrt(1.0-bb*bb)
      return self.mass*(gam-1) 
           
    def averageRadialVelocity(self):  # in cUnit
      bb=0.
      if (len(self.data) < 2): 
        return 0.
      for ee in self.data:
        bb+=numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3])
      bb/=len(self.data)
      return bb

    def averageRadius(self):  # in keV.. 
      rr=0.
      if (len(self.data) < 2): 
        return 0.
      for ee in self.data:
        rr+=numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
      rr/=len(self.data)
      return rr      

    # Propagation, by field integration, while is passing by.  
    def propagateWithBeam(self, phiFromBunch, staticBField, mXPipe, mYPipe, totalQ, tokenTraj ):
        # Check the arguments.  Expect a Real scalar field 
      argName1=str(phiFromBunch.__class__)
      if (argName1.find("Real_scalar_field") == -1): 
        print " Electron Flock: Wrong 1st argument type, expect a Real_scalar_field, type is  ", \
	    argName1
	print " Fatal Error "     
	sys.exit() 
	    
      argName2=str(staticBField.__class__)
      if (argName2.find("array") == -1): 
        print " Electron Flock: Wrong 2nd argument type, expect an array, type is  ", \
	    argName2
	print " Fatal Error "     
	sys.exit() 
     # all particles are treated separatly.
     # Instantiate an Runge Kutta Propagator
      
      myRK = ECloudPy.RKIntegrator(False);
      myRK.setDynamicRelativistic(True);
      myRK.setMaximumXBeamPipe(5.5e-2);
      myRK.setMaximumYBeamPipe(2.5e-2);
      myRK.setBFieldStaticCmp(staticBField[0], 0);  
      myRK.setBFieldStaticCmp(staticBField[1], 1);  
      myRK.setBFieldStaticCmp(staticBField[2], 2);
      # second argument is probably wrong if frequency scale is not set to 1. 
      myRK.setUnits(totalQ*0.1, 1.0); 
      dataTmp=[]
      dataBPTmp=[]
      nCnt=1
      myRK.setDebugOff()
      self.nBad=0
      for ee in self.data:
        tOff = -1.0*ee[4]/synergia.PH_MKS_c
	tFinal = numpy.array([0.],'d')
	rIn = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
        if (nCnt < 4):
	  fNameNum="PyTrajElecNum%d" % nCnt
	  fName=tokenTraj+fNameNum+".txt"
	  print "... declaring file name", fName 
	  myRK.reOpenTrajectoryFile(fName) 
	myRK.propThroughBunch(ee,phiFromBunch, tOff, tFinal)
	if (myRK.gotPropagationError()):
	  self.nBad=self.nBad+1
	  continue 
        if (nCnt < 4):
	  myRK.closeTrajectoryFile()
	rOut = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
      #  they all land at slightly different time..
	ee[6]=tFinal[0] 
        if (myRK.reachedBeamPipe()): 
#	  print " Reached BeamPipe!!!  ", rOut 
          dataBPTmp.append(ee)
        else:
	  dataTmp.append(ee)
#	  print " In vaccum Radius out for electron ", nCnt, " = ", rOut 
	nCnt=nCnt+1	
#	print " At electron  ", nCnt, " = ", rOut 
      # End loop on electrons. 	  
      self.data = dataTmp
      self.dataBP = dataBPTmp
     # Propagation, by field integration, in between bunches   
    def propagateNoBeam(self, staticBField, deltaT):
	    
      argName2=str(staticBField.__class__)
      if (argName1.find("array") == -1): 
        print " Electron Flock: Wrong 2nd argument type, expect an array, type is  ", \
	    argName2
	print " Fatal Error "     
	sys.exit() 
     
    def pyplotEk1(self, token, reachedBP, maxEk):
    
      import pylab
      nn=0
      fName=""
      if reachedBP:
        nn=len(self.dataBP)
	if (nn<2):return 
        eks=numpy.zeros(nn)
        ephis=numpy.zeros(nn)
        ii=0
        for ee in self.dataBP:
          beta=numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3] + ee[5]*ee[5])
          gam = 1./numpy.sqrt(1.0-beta*beta)
	  ek=self.mass*(gam-1)
	  if (ek > maxEk): ek=maxEk
          eks[ii]=ek
	  ephis[ii]=numpy.arctan2(ee[3], ee[1]) 
          ii=ii+1
        pylab.hist(eks, 100, log=True, bottom=1.)
	fName="EK_"+token+"BP.pdf"
        pylab.savefig(fName)
        pylab.clf()
        pylab.hist(ephis, 50, log=False, bottom=1.)
	fName="Phis_"+token+"BP.pdf"
        pylab.savefig(fName)
        pylab.clf()
	
      else:
        nn=len(self.data)
	if (nn<2):return 
        eks=numpy.zeros(nn)
        ephis=numpy.zeros(nn)
        ii=0
        for ee in self.data:
          beta=numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3] + ee[5]*ee[5])
          gam = 1./numpy.sqrt(1.0-beta*beta)
	  ek=self.mass*(gam-1)
	  if (ek > maxEk): ek=maxEk
          eks[ii]=ek 
	  ephis[ii]=numpy.arctan2(ee[3], ee[1]) 
          ii=ii+1
        pylab.hist(eks, 100, log=True, bottom=1.)
	fName="EK_"+token+"VC.pdf"
        pylab.savefig(fName)
        pylab.clf()
        pylab.hist(ephis, 50, log=False, bottom=1.)
	fName="Phis_"+token+"VC.pdf"
        pylab.savefig(fName)
        pylab.clf()
    
