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
      self.totalCharge=0. # in Coulomb, for electrons in this flock
      self.mxPipe=5.5e-2 # in meters. 
      self.myPipe=2.5e-2 # in meters.
      self.numTrajDebug=4 # Number of trajectory dumped on file for analysis
      # Variable from controlling multibunch simulation. 
      # a complete single bunch phase consists of 
      # (i) propagating the electon in vaccum through a bunch crossing
      # (ii) collecting secondary emission from the wall, during bunch crossing.
      # (iii) propagating those through the bunch   
      # (iv) progating of electrons between bunches. 
      # (v) collecting secondary emission from the wall, between bunch crossing.
      # (vi) propagating those 
      # (vii) generating electron from residual gas, 
      # (viii) progating these electrons through the bunch and between the bunch.
      # At the end of one phase, all electron clocks should be alinged at one bunch crossing,
      # all electron should have their clock at finalBunchClock, which is the time between 
      # bunches.   
      #  
      self.bunchSpacing=18.8e-9 # in seconds(default FNAL r.f. freq. of 53.3 MHz
      self.betaProtonBunch=0.994478 # velocity of the Proton Bunch 
      self.electronCount=0. # we need a counter to uniquely tag the electron. 
      # Some properties of the proton bunch are convenient to be kept in memory here
      self.numIonsFromGasPerLength=25000
      self.totalChargeProtonBunch=3e11 # in Units appropriate to compute e.m. force from potential 
      self.xOffsetProtonBunch=0.
      self.yOffsetProtonBunch=0.
      self.zOffsetProtonBunch=0.
      self.xWidthProtonBunch=0.005
      self.yWidthProtonBunch=0.0035
      self.protonBunchLength=0.45
      self.bunchNumber=0
      
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

    def setBunchSpacing(self, bb):
      self.bunchSpacing=bb
      
    def setBetaProtonBunch(self, bb):
      self.betaProtonBunch=bb
      
#   
    def setNumIonsFromGasPerLength(self, n):  
      self.numIonsFromGasPerLength=n
      
    def setXOffsetProtonBunch(self, xo):  
      self.xOffsetProtonBunch=xo
      
    def setYOffsetProtonBunch(self, yo):  
      self.yOffsetProtonBunch=yo
      
    def setZOffsetProtonBunch(self, zo):  
      self.zOffsetProtonBunch=zo
    
    def setXWidthProtonBunch(self, xw):  
      self.xWidthProtonBunch=xw
      
    def setYWidthProtonBunch(self, yw):
      self.yWidthProtonBunch=yw
     
    def setProtonBunchLength(self, bl): 
      self.protonBunchLength=bl
      
    def setTotalChargeProtonBunch(self, totalQ):
      self.totalChargeProtonBunch=totalQ
      
    def addElectron(self, anElectron):
    # No check of argument, my excuse, no time for this.. 
      self.data.append(anElectron)

    def clear(self):
      self.data = []
      self.dataBP = []
   
    def resetBunchNumber(self):
      self.bunchNumber=0
   
    # Adding electron.  over the length of the bunch. 
    def addFromGas(self, prescaleFact):
     
      gamProtons = numpy.sqrt(1.0/(1.0-self.betaProtonBunch*self.betaProtonBunch))
      massElecOMassProton = synergia.PH_NORM_me/synergia.PH_NORM_mp
      	   
# prescale factor of 100, for now...     
      num_electrons =  int(self.numIonsFromGasPerLength*prescaleFact)
    
      aGaussX=numpy.random.normal(self.xOffsetProtonBunch, self.xWidthProtonBunch, 100)
      aGaussY=numpy.random.normal(self.yOffsetProtonBunch, self.yWidthProtonBunch, 100)
      # Choose here quasi steady state solution, i.e., flat distribution..
#      if steadyState:
#	aGaussZ=numpy.random.rand(100)
#	for k in range(100):
#	  aGaussZ[k]=aGaussZ[k]*bunchLength
#      else: 
      aGaussZ=numpy.random.normal(0., self.protonBunchLength, 100)
      EmaxElec = 2.0* self.mass * 1.0e3 * \
                (self.betaProtonBunch * self.betaProtonBunch/(gamProtons*gamProtons))/ \
                (1.0 + 2.0* gamProtons * massElecOMassProton + \
		massElecOMassProton*massElecOMassProton) # in eV
      for i in range(num_electrons):
        ii = i%100  # by chunks of 100 ... Could go faster ? 
        if ((ii == 0) & (i > 1)):
          aGaussX=numpy.random.normal(self.xOffsetProtonBunch, self.xWidthProtonBunch, 100)
          aGaussY=numpy.random.normal(self.yOffsetProtonBunch, self.yWidthProtonBunch, 100)
#          if steadyState:
#	    aGaussZ=numpy.random.rand(100)
#	    for k in range(100):
#	      aGaussZ[k]=aGaussZ[k]*bunchLength
#	  else: 
          aGaussZ=numpy.random.normal(0., self.protonBunchLength, 100)
	  #Assume for now average energy of ~ 3 eV, Landau distributed.
        
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
	tClock=(aGaussZ[ii]+self.zOffsetProtonBunch)/synergia.PH_MKS_c
        thisElec = numpy.array([aGaussX[ii], betax, aGaussY[ii], betay, \
                              aGaussZ[ii], betaz, tClock, 0., self.electronCount])
	self.electronCount += 1.		      
#      print "an Elec", thisElec
#        if (numpy.abs(thisElec[2]) > 0.023):
#	  print " Large Gaussian tails in Y " 
#	  sys.exit()
	  # O.K., did not happened...  
        self.addElectron(thisElec)
    
      self.setTotalCharge(synergia.PH_MKS_e * self.numIonsFromGasPerLength)
    
    def numInVaccum(self):
      return len(self.data) 
      
    def numLeftToTrackBC(self): # Stands for left to track between Bunch Crossing... 
      nToTrack=0
      for ee in self.data:
        if (ee[6] < self.bunchSpacing): nToTrack=nToTrack+1
      return nToTrack 
      
    def numLeftToTrack(self, tMax):
      nToTrack=0
      for ee in self.data:
        if (ee[6] < tMax): nToTrack=nToTrack+1
      return nToTrack 
      
    def numReachedBeamPipe(self):
      return len(self.dataBP) 
    
    def electronCounter(self):
      return int(self.electronCount)
      
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
      
    def averageClockBP(self):  # Hitting the Beam Pipe... in ns
      bb=0.
      if (len(self.dataBP) < 2): 
        return 0.
      for ee in self.dataBP:
        bb+=ee[6]
      bb/=len(self.dataBP)
      return bb*1.0e9 # in ns  
      
    def averageTimeOfFlight(self):  # Hitting the Beam Pipe... in ns
      bb=0.
      if (len(self.data) < 2): 
        return 0.
      for ee in self.data:
        bb+=ee[7]
      bb/=len(self.data)
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
      # Before paying the price of instantiating a Integrator, look if we 
      # something totrack 
      if (self.numLeftToTrackBC()==0):
        return
      physSize= phiFromBunch.get_physical_size()
      deltaTMax=physSize[2]/(self.betaProtonBunch*synergia.PH_MKS_c) - 1.0e-11 
      # 10 ps safety margin...
#      print " deltaTMax in propagateWithBeam ", deltaTMax
      nToTrackHere=0
      for ee in self.data:
        if (ee[6] < deltaTMax): nToTrackHere=nToTrackHere+1
      
      if (nToTrackHere == 0): return
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
#       Not needed, as the propagation will occur through the end of the bunch.
      for ee in self.data:
        if (ee[6] > deltaTMax):
#	  print " Skipping electron ", int(ee[8]), " at time ", ee[6], " z = ",  ee[4]
	  dataTmp.append(ee)
	  continue
#	print " Tracking electron ", int(ee[8]), " at time " , ee[6], "  z= " , ee[4]
        eeZOrig=ee[4]
	tOff = ee[6]
	tOffInit = tOff
        if steadyState:
	  ee[4] = -physSize[2]/2.
	  tOff=0.
	tFinal = numpy.array([0.],'d')
	rIn = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
        if (nCnt < self.numTrajDebug):
	  fNameNum="TrBTrajElecNum%d" % int(ee[8])
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
	    fNameNum="TrBBadTrajElecNum%d" % int(eeinit[8])
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
#	print " State for electron ", int(ee[8]), " at time " , ee[6], "  z= " , ee[4]
	if steadyState: 
	  ee[4] = ee[4] + eeZOrig; 
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
     # tokenTraj : a string should there be a request to dump trajectories.    
    def propagateNoBeam(self, staticBField, tokenTraj):
	    
      argName1=str(staticBField.__class__)
      if (argName1.find("array") == -1): 
        print " Electron Flock: Wrong first argument type, expect an array, type is  ", \
	    argName2
	print " Fatal Error "     
	sys.exit() 
      if (self.numLeftToTrackBC()==0):
        return
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
	deltaT = self.bunchSpacing - tOff
	if (numpy.abs(deltaT) < 1.0e-12):
	  dataTmp.append[ee]
	  continue
	tFinal = numpy.array([0.],'d')
	rIn = numpy.sqrt(ee[0]*ee[0] + ee[2]*ee[2])
        if (nCnt < self.numTrajDebug):
	  fNameNum="NBTrajElecNum%d" % nCnt
	  fName=tokenTraj+fNameNum+".txt"
	  myRK.reOpenTrajectoryFile(fName)
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
#	  print " Reached BeamPipe in between bunch!!!, time   ", ee[6]
          dataBPTmp.append(ee)
        else:
	  ee[6]=self.bunchSpacing + 2.0e-12 # Done, add a pico or two for good measure.  
	  dataTmp.append(ee)
#	  print " In vaccum Radius out for electron ", nCnt, " = ", rOut 
	nCnt=nCnt+1	
#	print " At electron  ", nCnt, " = ", rOut 
      # End loop on electrons. 	  
      self.data = dataTmp
      self.dataBP = dataBPTmp
    
    def AddFromWall(self, seyRate, materialIndex):
     #  Take the electron stuck on the beam pipe 
     
       numIons=len(self.data)
       for ee in self.dataBP: # Looping on electron hitting the Beam Pipe
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
	 for k in range(nSecTmp):
	 # not sure at all about the geometry!... 
	   betaR=stuff[0][k] # First element in the tuple is the normal, i.e. radial for us
	   betaZ=numpy.sqrt(stuff[1][k]*stuff[1][k]+stuff[2][k]*stuff[2][k])
	  # Information loss here: the phi angle with repsect to the normal to 
	  # is missing.. Not sure this makes sense.. 
	   betaX=-1.0*betaR*numpy.cos(incPhiRadial) # away from the wall
	   betaY=-1.0*betaR*numpy.sin(incPhiRadial)
	   # at any rate, not quite correct for elliptical geometry. 
	   dtt=1.0e-15 # a symbolic one femto delay...
	   tClock = ee[6]+dtt	  
	   # Shift the elctron a bit such that he is no longer in the metal..
	   if (ee[0] < 0.): 
	   	ee[0] += 1.0e-4
	   else:
	        ee[0] -= 1.0e-4
		
	   if (ee[2] < 0.): 
	   	ee[2] += 1.0e-4
	   else:
	        ee[2] -= 1.0e-4
	   # Introducing a bug-feature here to test:   going back in time 
	   # ee[6] = 0.5e-9 
	   # O.K., now the secondary senses the bunch...	
           thisElec = numpy.array([ee[0], betaX, ee[2], betaY, \
                              ee[4], betaZ, tClock, dtt, self.electronCount])
	   self.electronCount += 1.		      
#	   print "from eIn ", incEnergy[0], \
#	         "Added Wall electrons betaX=", betaX, " betaY=", betaY, \
#	         " betaZ=", betaZ 		       
           self.addElectron(thisElec)
	   numIons=numIons+1
         # end loop secondary electrons
       # end loop primary electrons 
       self.setTotalCharge(synergia.PH_MKS_e * numIons)
       self.dataBP=[] # we got rid of those.. NO re-DIFFUSE ELECTRONS !? 
    
    
    def propagateOneCrossing(self, prescaleFactor, phiFromBunch, bField, tokenCase):
    
      self.resetClockForNextBunch()
      physSize= phiFromBunch.get_physical_size()
      deltaTMax=physSize[2]/(self.betaProtonBunch*synergia.PH_MKS_c) - 1.0e-11
      print " PropagateOneCrossing, deltaTMax ", deltaTMax
      # We propagate the electron left over from the previous bunch crossing, 
      # though the current bunch crossing, until we are running out of regenerated electron, 
      # for the bunch crossing time..
      nIterPhase1=0
      steadyState=True  
      while (self.numLeftToTrack(deltaTMax) > 0):
        self.propagateWithBeam(steadyState, phiFromBunch, bField, self.totalChargeProtonBunch, "")
	print " SteadyState BunchCrossing, iter ", nIterPhase1,\
	      " numElectrons ", len(self.data), " Reached Pipe " , len(self.dataBP)
	print " Average Clock on Beam Pipe ", self.averageClockBP()      
	self.AddFromWall(2.0, 1)
	steadyState=False # Late electron do not see the whole bunch...
	nIterPhase1+=1

      if (len(tokenCase) > 2):
        fName="ElectronFlockBC_" + tokenCase + "_Crossing_%d" % self.bunchNumber + ".txt" 
        self.writeToASCIIFile(fName)

      # Now add the newly created electron over the length of bunch crossing. 	   
      self.addFromGas(prescaleFactor)
      print "  Added from Gas interaction, numElectrons ", len(self.data)
      #
      # 
      nIterPhase2=0  
      while (self.numLeftToTrack(deltaTMax) > 0):
        self.propagateWithBeam(False, phiFromBunch, bField, self.totalChargeProtonBunch, "")  
	self.AddFromWall(2.0, 1)
	print "  After from Wall Sec Emission, iter ", nIterPhase2, " numElectrons ", len(self.data)
	nIterPhase2+=1
	
      # Now we drift between bunches.. Again, regenerating electrons.. 
      nIterPhase3=0  
      while (self.numLeftToTrackBC() > 0):
        self.propagateNoBeam(bField, "")
	self.AddFromWall(2.0, 1)
	print "  After from Sec emission, between bunches iter ", \
	       nIterPhase3, " numElectrons ", len(self.data)
	nIterPhase3+=1
	
      print " PropagateOneCrossing, Ending with ", len(self.data)
      if (len(tokenCase) > 2):
        fName="ElectronFlockDR_" + tokenCase + "_Crossing_%d" % self.bunchNumber + ".txt" 
        self.writeToASCIIFile(fName)
      self.bunchNumber+=1
      
    def resetClockForNextBunch(self):
      for ee in self.data:
	  ee[6]=0.
	  ee[7]=0.
    # Unbiased resampling 
    def reSample(self, fraction):
      dataTmp=[]
      rr=numpy.random.rand(len(self.data))
      i=0
      for ee in self.data:
          if (rr[i] < fraction): 
	    dataTmp.append(ee)
	  i+=1
      self.data=dataTmp 
   
    
    def writeToASCIIFile(self, fName):
      output=open(fName,'w')
      output.write("n x bx y by z bz t dt \n")
      for ee in self.data:
        line = " %d" % int(ee[8])
	for k in range(3):
	  line += " %10.5g " % (1000.*ee[2*k])  
	  line += " %10.5g " % ee[2*k+1]
	line += " %10.5g " % (1.0e9*ee[6])    
	line += " %10.5g \n " % (1.0e9*ee[7])
        output.write(line)
      output.close()
       
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
        
    
