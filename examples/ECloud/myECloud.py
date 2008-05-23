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

 
if ( __name__ == '__main__'):

    t0 = time.time()
    current = 0.5
    kinetic_energy = 8.0  # Injection energy.. 
    mass = synergia.PH_NORM_mp
    massElecOMassProton = synergia.PH_NORM_me/synergia.PH_NORM_mp
    massE = synergia.PH_NORM_me
    massEev = 1.0e9*synergia.PH_NORM_me
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6 # Angular frequency 
    part_per_cell = 1
    pipe_radius = 0.0254 # one inch beam pipe
    kicks_per_line = 10 # Number of iteration on space charge field kick 
    gridnum = 32 
    griddim = (gridnum,gridnum,gridnum)
    num_particles = griddim[0]*griddim[1]*griddim[2] * part_per_cell
    solver = "3D"
    
    print "num_particles =",num_particles
    print "We will use a", solver, "solver"
    
    xwidth=0.0032   # 5 pi mm mrad, beta = 55 m., gamma = 9.  
    ywidth=0.0056   # who know, emittance is a bit bigger.. 
    xpwidth=0.0049608  # not sure it matters for now.. 
    xOffset = 0.002
    yOffset = 0.005
    rx=0.85440 
    dpop = 1.0e-20
    bunchLength = 0.45*c_light*1.0e-9 # in meters 
    protonPerBunch = 2.0e12

    # generate electron. Get the Tech-X Cross section for proton on gas.. 
    
   # Extracted from Spentz prototype..

    # These are output variables
    # set velocity
    oneOGam = proton_mass/(kinetic_energy + proton_mass)
    betaProtons=(1.0 - oneOGam*oneOGam)**0.5
    gamProtons = 1.0/oneOGam
    print " beta protons " , betaProtons
    vProtons = betaProtons*c_light
   
#    ee = synergia.Error_eater()
#    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.mad"),"channel",kinetic_energy,
                        scaling_frequency)
			
    gourmet.insert_space_charge_markers(kicks_per_line)

    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)
    print "Mass of particle", mass
    print "Beam Beta", beam_parameters.get_beta()
    print "Beam Gamma", beam_parameters.get_gamma()
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 
    print "Betagamma and inverse betagamma",betagamma,1./betagamma
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = -rx)
    beam_parameters.y_params(sigma = xwidth, lam = xpwidth * pz,r = rx)
    sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
    beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)

    sys.stdout.flush()
    
    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)
    
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
    print "elapsed time =",time.time() - t0

#
# How do I get the sigma of the bunch after propagation through the lattice 
# put the injected parameters for now.

    means,sigs = synergia.get_spatial_means_stds(bunch)
    sigXNow = sigs[0]
    sigYNow = sigs[1]
    bunchLength = sigs[2] 
    print " Bunch sigmas, x =  ", sigXNow, " y = ", sigYNow, " z ", bunchLength

    this_vel = numpy.array([vProtons])
    flagIon=1
    flagGas=5 # O2, ignore the hydrogen in the water..  
      
    myCross = txphysics.txionpack.get_sigma_impact_array(this_vel, flagGas, flagIon)
    # This call fills the variables numColls and which_gs 
    xsec_Mbarn=myCross/m_sqr_to_Mbarn
    vacuumPressure=2.0e-7
    vacuumTemperature=305
    ion_elecs = protonPerBunch*3.28 * xsec_Mbarn * vacuumPressure * 294 / vacuumTemperature
    # factor 3.28 = torr to Pascal (133.32) * Avogadro * / (R=8.314 J. mol^-1.K^-1 * 294 K) 
    # Assume the cross section are for molecule, not atomic.. 
    # This quantity is the linear density, electrons and ions per meter
    print "Ionization production: ",ion_elecs, " from ", protonPerBunch, " protons"
    print "Production cross section [Mbarns]:", xsec_Mbarn
         
#   Create the flock of electron.
#
    mEl = ElectronFlock()
#    
# Prescale factor for creating to the flock.. 
#
    prescaleFact = 0.01; 
        
    mEl.addFromGas(ion_elecs, prescaleFact, betaProtons, \
                    xOffset, yOffset, sigXNow, sigYNow, bunchLength)
    
    print " Number of electron created ", mEl.numInVaccum()
    print " Average kinetic energy for my electron flock ", mEl.averageKineticEnergy()
    print " Average radius for this flock ", mEl.averageRadius()
    print " Average radial velocity ", mEl.averageRadialVelocity()
    
# How do I get the potential from the bunch after propagation?
# define a new scalar and re-populate with charge from the bunch ?  
# repopulate with charge? This looks like a waste of time !. 
# get the sigma of the bunch.  For now, assume it did not changed by the 
# by the propagation...      
    sizeNow = (5.0*sigXNow, 5.0*sigYNow, 3.0*bunchLength)
# set the offset.. ?? Keep it centered     
    rhoAfterProp = Real_scalar_field(griddim,sizeNow,(0.0,0.0,0.0))
# drop the bunch charge onto the grid..    
    total_charge = deposit_charge_cic(rhoAfterProp,bunch.get_store(),0)
# Recompute the potential
#    griddimA = numpy.array([griddim[0], griddim[1], griddim[2]])
    fftwh = Fftw_helper(griddim, False)
    phiAfterProp = solver_fftw_open(rhoAfterProp,fftwh,0)
#
# Now propagate the electron through this potential, with a static Bfield 
# of 1 kG oriented along the Y axis..     
#
    bField = numpy.array([0., 0.1, 0.], 'd')
    mEl.propagateWithBeam(phiAfterProp, bField, mXPipe=5.5e-2, mYPipe=2.5e-2)
    print " Number of electron in Vacuum ", mEl.numInVaccum()
    print " Number of electron Reached Beam Pipe ", mEl.numReachedBeamPipe()
    print " Average kinetic energy for my electron in vaccum ", mEl.averageKineticEnergy()
    print " Average radius... ", mEl.averageRadius()
    print " Average radial velocity ", mEl.averageRadialVelocity()
    
# End of electron loop.

#    dimpact = synergia.Diagnostics_impact_orig("channel_impact_open")
#    d0 = synergia.Diagnostics_impact_orig("channel0current")
    print " One bunch propagated though the potential of the bunch, enough for today.. " 
    sys.exit()

    import pylab

#    pylab.plot(d0.s,d0.std[:,synergia.x],'gx',label='no space charge')
#    pylab.xlabel('s (m)')
#    pylab.ylabel('std<x> (m)')
#
#    pylab.plot(dimpact.s,dimpact.std[:,synergia.x],'o',label='impact')
#    pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],'r+',markersize=15.0,label='fish')
#    pylab.legend(loc=0)
#    pylab.show()
