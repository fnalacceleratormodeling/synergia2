#!/usr/bin/env bwpython

import numpy
import time
import math
import os
import sys

import synergia
import s2_fish
from s2_solver_fftw import *
from s2_containers import *
from s2_deposit import *
#
# some global constants

c_light = 299792458.0
# in (m/s)
electron_mass_kg = 9.10938188e-31
# in (kg)
proton_charge = 1.60217646e-19
# in C
epsilon_0 = 8.85418782e-12
# electric constant [m^{-3} kg^{-1} s^4 A^2] or [F/m]
proton_mass = 0.9382723
electron_mass = 0.51099906e-3
# in GeV
m_sqr_to_Mbarn = 10**-22
mpoverme = 1836
# proton mass over electron mass


from mpi4py import MPI
import txphysics.txionpack
import ECloudPy
from electronFlock import ElectronFlock
from testMPISimple import testMPISimple
#import plotPotential
 
if ( __name__ == '__main__'):

#    fNameLog="LogFileECloud1NodeV2_%d" % MPI.COMM_WORLD.Get_rank() + ".lis"
#    outputLog=open(fNameLog,'w')

    print " My MPI.COMM_WORLD.Get_rank() is ", MPI.COMM_WORLD.Get_rank()
#    outputLog.write(" My MPI.COMM_WORLD.Get_rank() is %d \n" % MPI.COMM_WORLD.Get_rank())
#    outputLog.flush()
#    outputLog.close()
#    sys.exit() 
#    print " My Processor Name is ", MPI.Get_processor_name()
#    outputLog.write(" My Processor Name is " + MPI.Get_processor_name() + " \n")
# Setting the numpy random seeds.. 
# Start different seeds based on the rank.  Crummy!!! Not 
    junk1=numpy.random.rand(1)
    numThrow = int(100*junk1[0])*MPI.COMM_WORLD.Get_rank() + 1
#    outputLog.write(" Number of throw %d \n" % numThrow)
    for k in range(numThrow):
      aa=numpy.random.rand(100)
    
        
#    mT=testMPISimple()
#    if (MPI.COMM_WORLD.Get_rank() == 0):
#      print " Generating .... " 
#    mT.generate()
#    if (MPI.COMM_WORLD.Get_rank() == 0):
#      print " Analysis .... " 
#   mT.analyze()   
#    if (MPI.COMM_WORLD.Get_rank() == 0):
#      print " Collection .... " 
#    mT.sumOut()
#    mT.spitOut()
#    mT.spitOutFast()
#    MPI.Finalize()
#   sys.exit()
       
    t0 = time.time()
    current = 0.5
    kinetic_energy = 8.0  # Injection energy.. 
    mass = synergia.PH_NORM_mp
    massElecOMassProton = synergia.PH_NORM_me/synergia.PH_NORM_mp
    massE = synergia.PH_NORM_me
    massEev = 1.0e9*synergia.PH_NORM_me
    charge = 1.0
    initial_phase = 0.0
#    scaling_frequency = 10221.05558e6 # Angular frequency 
    scaling_frequency = 2.0e8 # Angular frequency 
    scaling_frequency = synergia.PH_MKS_c/(2.0*math.pi) # Angular frequency 
    part_per_cell = 1
    pipe_radius = 0.0254 # one inch beam pipe, Radius!
    kicks_per_line = 10 # Number of iteration on space charge field kick 
    gridnum = 64 # 64 is the default... 	 
    griddim = (gridnum,gridnum,2*gridnum)
    num_particles = griddim[0]*griddim[1]*griddim[2] * part_per_cell
    solver = "3D"
    
    if MPI.COMM_WORLD.Get_rank() == 0:
      print "num_particles =",num_particles
      print "We will use a", solver, "solver"
    
    xwidth=0.0050   # 25 pi mm mrad, beta = 55 m., gamma = 9.  
    ywidth=0.0055   # who know, emittance is a bit bigger.. 
    xpwidth=0.0049608  # not sure it matters for now.. 
    ypwidth=0.00768  # random crap for now, matching conditions to be studied later
    xOffset = 0.000
    yOffset = 0.0000001
    rx=0.85440 
    dpop = 1.0e-20
    bunchLength = 1.5*c_light*1.0e-9 # in meters 
    protonPerBunch = 3.1e11 # Nominal is 3.1 10^11 particle per bunch 
    tokenCase="V3e11ppbMediumPNegLW16" # MediumP means 2.0e-8 
    
    # generate electron. Get the Tech-X Cross section for proton on gas.. 
    
   # Extracted from Spentz prototype..

    # These are output variables
    # set velocity
    oneOGam = proton_mass/(kinetic_energy + proton_mass)
    betaProtons=(1.0 - oneOGam*oneOGam)**0.5
    gamProtons = 1.0/oneOGam
    if MPI.COMM_WORLD.Get_rank() == 0:
      print " beta protons " , betaProtons
    vProtons = betaProtons*c_light
   
#    ee = synergia.Error_eater()
#    ee.start()
# Should I do this on every node ? 
# Guess so... 
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.mad"),"channel",kinetic_energy,
                        scaling_frequency)
			
    gourmet.insert_space_charge_markers(kicks_per_line)

    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=0)
# note: last argument, we want 6D Gaussian beam... 					 
    if MPI.COMM_WORLD.Get_rank() == 0:
      print "Mass of particle", mass
      print "Beam Beta", beam_parameters.get_beta()
      print "Beam Gamma", beam_parameters.get_gamma()
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 
    if MPI.COMM_WORLD.Get_rank() == 0:
      print "Betagamma and inverse betagamma",betagamma,1./betagamma
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = -rx)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = rx)
    sigma_z_meters = bunchLength
    if MPI.COMM_WORLD.Get_rank() == 0:
      print " sigma_z_meters ", sigma_z_meters
    beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)
#    print " And Quit " 
#    sys.exit()
    
    sys.stdout.flush()
    
    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)
    means,sigs = synergia.get_spatial_means_stds(bunch)
    sigXNow = sigs[0]
    sigYNow = sigs[1]
    bunchLength = sigs[2] 
    if MPI.COMM_WORLD.Get_rank() == 0:
      print " Bunch sigmas, before propagation, x =  ", \
            sigXNow, " y = ", sigYNow, " z ", bunchLength
    
    line_length = gourmet.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    diag = synergia.Diagnostics(gourmet.get_initial_u())
    kick_time = 0.0
    
    if solver == "3D" or solver == "3d":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_s2_fish=True,periodic=True)
    elif solver =="2D" or solver == "2d":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_gauss=True)
#    print "elapsed time =",time.time() - t0
    unitsBunch = bunch.get_store().get_units()
#    print " Units X or Y for the bunch ", unitsBunch[0]  
#    print " Units Z for the bunch ", unitsBunch[4]    
#    print " And quit "
#    sys.exit()  
#
# How do I get the sigma of the bunch after propagation through the lattice 
# put the injected parameters for now.

    means,sigs = synergia.get_spatial_means_stds(bunch)
    sigXNow = sigs[0]
    sigYNow = sigs[1]
    bunchLength = sigs[2]
    print " Bunch sigmas, Rank ", MPI.COMM_WORLD.Get_rank(), ", sigYNow, y = ", sigYNow, " z ", bunchLength
# Take this bunch and compute a potential from it.     
# set the offset.. ?? Keep it centered     
#
    maxBeanPipeInX=0.055
    maxBeanPipeInY=0.025 # Radius..
    sizeNow = (2.0*maxBeanPipeInX+0.0005, 2.0*maxBeanPipeInY+0.0005, 6.0*bunchLength)
    rhoAfterProp = Real_scalar_field(griddim,sizeNow,(0.0,0.0,0.0))
# drop the bunch charge onto the grid..    
    total_charge = deposit_charge_cic(rhoAfterProp,bunch.get_store(),0)
    print " ..Total Charge Rank ", MPI.COMM_WORLD.Get_rank(), " Q = ", total_charge 
# Recompute the potential
#    griddimA = numpy.array([griddim[0], griddim[1], griddim[2]])
    fftwh = Fftw_helper(griddim, False)
#    print " ..Got FFTW Helper, Rank   ", MPI.COMM_WORLD.Get_rank()
    phiStripedAfterProp = solver_fftw_open(rhoAfterProp,fftwh, 0, True)
    phiAfterProp = Real_scalar_field(griddim,sizeNow,(0.0,0.0,0.0))
#    print " ..Did Solve , Rank   ", MPI.COMM_WORLD.Get_rank()
    i_lower = phiStripedAfterProp.get_points().get_dim0_lower()
    i_upper = phiStripedAfterProp.get_points().get_dim0_upper()
#    print " Rank ", MPI.COMM_WORLD.Get_rank(), " From myECloud,.. dim lower/upper ", i_lower , " / " , i_upper
     
    broadcast_Phi(phiStripedAfterProp, phiAfterProp, i_lower, i_upper)
#    print " ..Did broadcast    ", MPI.COMM_WORLD.Get_rank()
    
# put a barrier, as we want to test the potential NOW.     
#    MPI.COMM_WORLD.Barrier()
#    print " ..Did Barrier for the hell of it    ", MPI.COMM_WORLD.Get_rank()
    # Diagnostics for parallel processing: dump potential values from two different nodes. 
    #
    testPotentialMPI=False
    if (testPotentialMPI):
      print " On Node ", MPI.COMM_WORLD.Get_rank(), " getting  physSizeDiag"
      physSizeDiag=phiAfterProp.get_physical_size(); 
      print " On Node ", MPI.COMM_WORLD.Get_rank(), " Defining loc"
      loc = numpy.array([physSizeDiag[0]/5., physSizeDiag[1]/3., physSizeDiag[2]/4.],'d')
      print " On Node ", MPI.COMM_WORLD.Get_rank(), " getting vv ", loc
      if (MPI.COMM_WORLD.Get_rank() < 10):
        vv = phiAfterProp.get_val(loc)/total_charge # UnNormalize to the total charge (units?) 
        print " On Node ", MPI.COMM_WORLD.Get_rank(), " getting derivatives ", loc
        vvDerX = phiAfterProp.get_deriv(loc,0)/total_charge # UnNormalize to the total charge (units?) 
        print " On Node ", MPI.COMM_WORLD.Get_rank(), " At location ", loc
        print " ....MPIRank", MPI.COMM_WORLD.Get_rank(), "....Potential value ", vv, " derivative, X ", vvDerX
      if (MPI.COMM_WORLD.Get_rank() == 0):
        fNamePhiTmp= "Potential_Grid_"
	ffN1="%d" % gridnum
	ffN2="%d" % MPI.COMM_WORLD.Get_size()
	fNamePhiTmp+= ffN1+"_Nodes_"+ffN2+".dat"
	phiAfterProp.write_to_file(fNamePhiTmp)
#
        
      MPI.Finalize()
      sys.exit()

    this_vel = numpy.array([vProtons])
    flagIon=1
    flagGas=5 # O2, ignore the hydrogen in the water..  
      
    myCross = txphysics.txionpack.get_sigma_impact_array(this_vel, flagGas, flagIon)
    # This call fills the variables numColls and which_gs 
    xsec_Mbarn=myCross/m_sqr_to_Mbarn
    vacuumPressure=2.0e-8
    vacuumTemperature=305
    ion_elecs = protonPerBunch*3.28 * xsec_Mbarn * vacuumPressure * 294 / vacuumTemperature
    # factor 3.28 = torr to Pascal (133.32) * Avogadro * / (R=8.314 J. mol^-1.K^-1 * 294 K) 
    # Assume the cross section are for molecule, not atomic.. 
    # This quantity is the linear density, electrons and ions per meter
    if (MPI.COMM_WORLD.Get_rank() == 0):  
      print "Ionization production: ",ion_elecs, " from ", protonPerBunch, " protons"
      print "Production cross section [Mbarns]:", xsec_Mbarn
         
#    
# Prescale factor for creating to the flock.. 
#
    prescaleFact = 1.0 
#   Create the flock of electron.
#
    mEl = ElectronFlock()
    mEl.setMaxXDimPipe(maxBeanPipeInX)
    mEl.setMaxYDimPipe(maxBeanPipeInY)

      
    totalQ = protonPerBunch*proton_charge/total_charge
#    plotPotential.plotPotentialX(phiAfterProp, totalQ)
#    ffNamePotData="Potential"+tokenCase+".dat"
#    phiAfterProp.write_to_file(ffNamePotData)
#    print " And quit after graphing the potential "
#    sys.exit()
# Now deal with electrons.. 
# Prepare for a run over a few bunches... 
# 
# Thus:
    physSize= phiAfterProp.get_physical_size()
# Create electrons due beam gas interactions.
    mEl.setNumIonsFromGasPerLength(ion_elecs)
    mEl.setBetaProtonBunch(betaProtons)
    mEl.setXOffsetProtonBunch(xOffset)     
    mEl.setYOffsetProtonBunch(yOffset) 
    mEl.setXWidthProtonBunch(sigXNow) 
    mEl.setXWidthProtonBunch(sigYNow) 
    mEl.setZOffsetProtonBunch(physSize[2]/2.) 
    mEl.setProtonBunchLength(bunchLength)
    mEl.setTotalChargeProtonBunch(totalQ)  
    mEl.setNumberTrajectoryDump(0)
    
    # Static Bfield
#    bYDip=0.091*kinetic_energy/8.0
#    bYDip = 0.5e-4  # Earth magnetic field. 
#    bYDip = 0.091  # at 8 Gev, Check..
    bYDip = 0.228  

    bField = numpy.array([0., bYDip, 0.], 'd')
    mEl.setUniformBField(bField)    
#
# Test the fancy magent models.. 4 means DIPQUADEDGES; 3 QUADRUPOLEEDGE
# See ECloudCC/MIDipoleEdge.h
#     
    mEl.setMagnetModel(0, bYDip) # Drift.. 

    if (MPI.COMM_WORLD.Get_rank()==0):
      print " Static By = " , bField[1] 
    
    numCrossing=25
    prescaleFactorForGas=0.2
    
    maxNumElec=100000
    fractReSampled=1.0
    if (MPI.COMM_WORLD.Get_rank() == 0): 
      print " Prescale factor for Gas = ",  prescaleFactorForGas
      print " Maximum number of electron after one bunch .. ", maxNumElec
    mEl.resetBunchNumber()  
    for kBunch in range(numCrossing):
#      print " Number of electrons Begining Bunch Crossing ", kBunch, " is  ", mEl.numInVaccum() 
      mEl.propagateOneCrossing(prescaleFactorForGas, phiAfterProp, tokenCase) 
      if (MPI.COMM_WORLD.Get_rank()==0):
        print " Total number of electrons at end Bunch Crossing ", kBunch, " is  ", mEl.numInVaccum()
        if (mEl.numInVaccum() > maxNumElec):
          fractReSampled=float(maxNumElec)/mEl.numInVaccum()
          print " Node ", MPI.COMM_WORLD.Get_rank(), " kBunch = ", kBunch,  " Resampling, fraction  ", fractReSampled
      
      fractReSampled=MPI.COMM_WORLD.Bcast(fractReSampled, root=0)  
      if (mEl.numInVaccum() > maxNumElec):
        mEl.reSample(fractReSampled)
      # 
      prescaleFactorForGas *= fractReSampled
      if (MPI.COMM_WORLD.Get_rank()==0):
        print " Current prescale for Gas for bunch ", kBunch, " is  ", prescaleFactorForGas
      # And this is weight ..
    MPI.Finalize()
    sys.exit()


   
