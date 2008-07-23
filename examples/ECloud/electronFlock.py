#!/usr/bin/env python
# Electron flock, 
import numpy
import synergia
import s2_fish
import ECloudPy
import txphysics.txegenelec
import sys

class ElectronFlock:
#
    def __init__(self):
    
      self.data = [] # intented to be a list of numpy arrays, 8 doubles.
                     # 6 for phase space (3*(x,beta_x), time, time of flight )
      self.dataBP = [] # Electrons stuck on the beam pipe,8 doubles.
      self.nBad=0  # Number of electron lost in inf. loop or other propagation errors. 
      self.mass = 1.0e6*synergia.PH_NORM_me # Energy units here are keV..
      self.useRK = 1
      self.timeSliceWb=0.2 # don't recall what that is about..
      self.totalCharge=0. # in Coulomb
      self.mxPipe=5.5e-2 # in meters. 
      self.myPipe=2.5e-2 # in meters.
      self.numTrajDebug=4 # Number of trajectory dumped on file for analysis 
#      
#   Initialize the GSL random number for Laundau distribution.
#
      ECloudPy.initGslRan()

    def setMaxXDimPipe(self, dMaxPipeX):
      self.mxPipe=dMaxPipeX

    def setMaxYDimPipe(self, dMaxPipeY):
      self.myPipe=dMaxPipeY

    def setTotalCharge(self, charge):
      self.totalCharge=charge
     
    def setNumberTrajectoryDump(self, n):
      self.numTrajDebug=n
      
    def addElectron(self, anElectron):
    # No check of argument, my excuse, no time for this.. 
      self.data.append(anElectron)

    def clear(self):
      self.data = []
      self.dataBP = []
   
    # Adding electron.  Over one meter.. 
    def addFromGas(self, numIons, prescaleFact, betaProtons, \
                   xOffset,  yOffset, xWidth, yWidth, bunchLength, steadyState): 
      gamProtons = numpy.sqrt(1.0/(1.0-betaProtons*betaProtons))
      massElecOMassProton = synergia.PH_NORM_me/synergia.PH_NORM_mp
      	   
# prescale factor of 100, for now...     
      num_electrons =  int(numIons*prescaleFact)
    
      aGaussX=numpy.random.normal(xOffset, xWidth, 100)
      aGaussY=numpy.random.normal(yOffset, yWidth, 100)
      # Choose here quasi steady state solution, i.e., flat distribution..
      if steadyState:
	aGaussZ=numpy.random.rand(100)
	for k in range(100):
	  aGaussZ[k]=aGaussZ[k]*bunchLength - 0.5*bunchLength
      else: 
        aGaussZ=numpy.random.normal(0., bunchLength, 100)
      EmaxElec = 2.0* self.mass * 1.0e3 * (betaProtons * betaProtons/(gamProtons*gamProtons))/ \
                (1.0 + 2.0* gamProtons * massElecOMassProton + \
		massElecOMassProton*massElecOMassProton) # in eV
      for i in range(num_electrons):
        ii = i%100  # by chunks of 100 ... Could go faster ? 
        if ((ii == 0) & (i > 1)):
      	  aGaussX=numpy.random.normal(xOffset, xWidth, 100)
          aGaussY=numpy.random.normal(yOffset, yWidth, 100)
          if steadyState:
	    aGaussZ=numpy.random.rand(100)
	    for k in range(100):
	      aGaussZ[k]=aGaussZ[k]*bunchLength - 0.5*bunchLength
	  else: 
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
                              aGaussZ[ii], betaz, 0., 0.])
#      print "an Elec", thisElec
#        if (numpy.abs(thisElec[2]) > 0.023):
#	  print " Large Gaussian tails in Y " 
#	  sys.exit()
	  # O.K., did not happened...  
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
           
    def averageKineticEnergyBP(self):  # Hitting the Beam Pipe... in keV.. 
      bb=0.
      if (len(self.dataBP) < 2): 
        return 0.
      for ee in self.dataBP:
        bb+=numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3] + ee[5]*ee[5])
      bb/=len(self.dataBP)
      gam = 1./numpy.sqrt(1.0-bb*bb)
      return self.mass*(gam-1) 
      
    def averageTimeOfFlightBP(self):  # Hitting the Beam Pipe... in ns
      bb=0.
      if (len(self.dataBP) < 2): 
        return 0.
      for ee in self.dataBP:
        bb+=ee[7]
      bb/=len(self.dataBP)
      return bb*1.0e9 # in ns  
      
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

    # Propagation, by field integration, while the bunch is passing by.
    # steadyState: if true, all electrons will be propagate through the entire length of the bunch
    #              The duration of the crossing is given by the physical length of the grid
    #              on which the potential is defined. 
    #              If false, the electrons are assumed to have been created from residual 
    #              gas interaction, and will only sense the remaining part of the bunch located
    #              behind this interaction point
    # phiFromBunch: The scalar array that described the static electric potnetial from the 
    #                 proton bunch estimated in the reference frame of the bunch.
    # staticBField: vector of dimension 3, the static magnetic field.
    # TotalQ:  the total Charge of the proton bunch, used to appropriately scale the 
    #            the potential
    # tokenTraj:  A string, used to compose a file for trajectories. 
    #
    def propagateWithBeam(self, steadyState, phiFromBunch, staticBField, totalQ, tokenTraj ):
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
      myRK.setMaximumXBeamPipe(self.mxPipe);
      myRK.setMaximumYBeamPipe(self.myPipe);
      myRK.setStepRatio(10.);
#      myRK.setToPositron(); # just for kicks, try positron instead of electrons
      myRK.setBFieldStaticCmp(staticBField[0], 0);  
      myRK.setBFieldStaticCmp(staticBField[1], 1);  
      myRK.setBFieldStaticCmp(staticBField[2], 2);
      myRK.setUnits(totalQ, 1.0); 
      dataTmp=[]
      dataBPTmp=[]
      nCnt=1
      myRK.setDebugOff()
      self.nBad=0
      # if in steady state, go through the enire bunch...
      # to be implemented by shifting the electron in front of the potential
      # or setting tOff...
      tOff = 0.;
      physSize= phiFromBunch.get_physical_size()
#       Not needed, as the propagation will occur through the end of the bunch.
#      deltaTMax=numpy.abs(zMostLeft[2])/synergia.PH_MKS_c
      for ee in self.data:
        eeZOrig=ee[4]
        if steadyState:
	  ee[4] = -physSize[2]/2.
	tOff = ee[4]/synergia.PH_MKS_c
	tOffInit = tOff
	tFinal = numpy.array([0.],'d')
	rIn = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
        if (nCnt < self.numTrajDebug):
	  fNameNum="TrBTrajElecNum%d" % nCnt
	  fName=tokenTraj+fNameNum+".txt"
	  myRK.reOpenTrajectoryFile(fName)
	eeInit=numpy.array(ee)
	myRK.propThroughBunch(ee,phiFromBunch, tOff, tFinal)
#	if (numpy.abs(myRK.getMaximumYBeamPipe() - 2.5e-2) > 1.0e-6):
#	  print " Bad boundary.... Quit !!! "
#	  sys.exit()  
        if (nCnt < self.numTrajDebug):
	  myRK.closeTrajectoryFile()
	if (myRK.gotPropagationError()):
	  self.nBad=self.nBad+1
	  print "..Bad Trajectory Y ..", eeInit[2]
	  if (self.nBad<self.numTrajDebug): 
	    fNameNum="TrBBadTrajElecNum%d" % nCnt
	    fName=tokenTraj+fNameNum+".txt"
	    myRK.reOpenTrajectoryFile(fName)
	    tFinal[0]=0.
	    myRK.propThroughBunch(eeInit,phiFromBunch, tOffInit, tFinal)
	    myRK.closeTrajectoryFile()
	  continue 
	rOut = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
      #  they all land at slightly different time..
	ee[6]=tFinal[0]
	ee[7]=tFinal[0]-tOff
	if steadyState: 
	  ee[4] = ee[4] + eeZOrig + physSize[2]/2.; 
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
     # Arguments
     # static Bfield : vector of dimension 3, the static magnetic field 
     # tFinalReq : estimate final time at the end of the integration period. 
     #          i.e., all electron will have their time coordinate near or at tFinal. 
     # tokenTraj : a string should there be a request to dump trajectories.    
    def propagateNoBeam(self, staticBField, tFinalReq, tokenTraj):
	    
      argName1=str(staticBField.__class__)
      if (argName1.find("array") == -1): 
        print " Electron Flock: Wrong first argument type, expect an array, type is  ", \
	    argName2
	print " Fatal Error "     
	sys.exit() 
      myRK = ECloudPy.RKIntegrator(False);
      myRK.setDynamicRelativistic(True);
      myRK.setMaximumXBeamPipe(self.mxPipe);
      myRK.setMaximumYBeamPipe(self.myPipe);
      myRK.setBFieldStaticCmp(staticBField[0], 0);  
      myRK.setBFieldStaticCmp(staticBField[1], 1);  
      myRK.setBFieldStaticCmp(staticBField[2], 2);
      myRK.setDebugOff()
      nCnt=1
      dataTmp=[]
      dataBPTmp=[]
      for ee in self.data:
#        tOff = -1.0*ee[4]/synergia.PH_MKS_c
        tOff = ee[6]
	tOffInit = tOff
	tFinal = numpy.array([0.],'d')
	rIn = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
        if (nCnt < self.numTrajDebug):
	  fNameNum="NBTrajElecNum%d" % nCnt
	  fName=tokenTraj+fNameNum+".txt"
	  myRK.reOpenTrajectoryFile(fName)
	eeInit=numpy.array(ee)
	deltaT=tFinalReq-tOff
	myRK.propBetweenBunches(ee, tOff, deltaT, tFinal)
#	if (numpy.abs(myRK.getMaximumYBeamPipe() - 2.5e-2) > 1.0e-6):
#	  print " Bad boundary.... Quit !!! "
#	  sys.exit()  
        if (nCnt < self.numTrajDebug):
	  myRK.closeTrajectoryFile()
	if (myRK.gotPropagationError()):
	  self.nBad=self.nBad+1
	rOut = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
      #  they all land at slightly different time..
	ee[6]=tFinal[0]
	ee[7]=tFinal[0]-tOff 
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
    
    def AddFromWall(self, seyRate, materialIndex):
     #  Take the electron stuck on the beam pipe 
     
       for ee in self.dataBP:
         incEnergy=numpy.array([0.], 'd')
	 betaInc = numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3] + ee[5]*ee[5])
	 gamInc = numpy.sqrt(1.0/(1.0-betaInc*betaInc))
	 incEnergy[0]=(gamInc-1.0)*synergia.PH_NORM_me*1.0e9 # In eV now 
	 betaR = numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3])
         incCosAngle=numpy.array([0.], 'd')
         incCosAngle[0]=betaR/betaInc
	 stuff=txphysics.txegenelec.nsec_array(incEnergy, incCosAngle, materialIndex)
	 # now unwrap what I got .. 
	 nSecTmp=stuff[0].size
	 incPhiRadial=numpy.arctan2(ee[3], ee[1])
	 numIons=len(self.data)
	 for k in range(nSecTmp):
	 # not sure at all about the geometry!... 
	   betaR=stuff[0][k] # First element in the tuple is the normal, i.e. radial for us
	   betaZ=numpy.sqrt(stuff[1][k]*stuff[1][k]+stuff[2][k]*stuff[2][k])
	  # Information loss here: the phi angle with repsect to the normal to 
	  # is missing.. Not sure this makes sense.. 
	   betaX=betaR*numpy.cos(incPhiRadial)
	   betaY=betaR*numpy.sin(incPhiRadial)
	   # at any rate, not quite correct for elliptical geometry. 
	   dtt=1.0e-6 # a symbolic one femto delay... 
           thisElec = numpy.array([ee[0], betaX, ee[2], betaY, \
                              ee[4], betaZ, ee[6]+dtt, dtt])
#	   print "from eIn ", incEnergy[0], \
#	         "Added Wall electrons betaX=", betaX, " betaY=", betaY, \
#	         " betaZ=", betaZ 		       
           self.addElectron(thisElec)
	   numIons=numIons+1
         # end loop secondary electrons
       # end loop primary electrons 
       self.setTotalCharge(synergia.PH_MKS_e * numIons)
       self.dataBP=[] # we got rid of those.. NO re-DIFFUSE ELECTRONS !? 
       print " Done added from Wall done, and quit  "
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
        etof=numpy.zeros(nn)
        ii=0
        for ee in self.dataBP:
          beta=numpy.sqrt(ee[1]*ee[1] + ee[3]*ee[3] + ee[5]*ee[5])
          gam = 1./numpy.sqrt(1.0-beta*beta)
	  ek=self.mass*(gam-1)
	  if (ek > maxEk): ek=maxEk
          eks[ii]=ek
	  ephis[ii]=numpy.arctan2(ee[3], ee[1])
	  etof[ii]=ee[7]*1.0e9
          ii=ii+1
        pylab.hist(eks, 100, log=True, bottom=1.)
	fName="EK_"+token+"BP.pdf"
        pylab.savefig(fName)
        pylab.clf()
        pylab.hist(ephis, 50, log=False, bottom=1.)
	fName="Phis_"+token+"BP.pdf"
        pylab.savefig(fName)
        pylab.clf()
        pylab.hist(etof, 100, log=True, bottom=1.)
	fName="Toof_"+token+"BP.pdf"
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
        
    
