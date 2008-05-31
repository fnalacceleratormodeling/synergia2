#!/usr/bin/env bwpython
#
# This simple example run the 3D solver for a spherical charge distribution, Gaussian 
# generated. It is shown that the resulting potential falls of like 1/r, as the grid has been 
# chosen to be 20 times the sigma of this 3D distribution. 
#  
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

from mpi4py import MPI
 
if ( __name__ == '__main__'):

    kinetic_energy = 8.0  # Injection energy..Irrelevant, it seems..  
    mass = synergia.PH_NORM_mp # mass of the proton, the bunch contains proton. 
    # Set the charge of the bunch.  Unspecified units 
    charge = 2.0  # unspecified units 
    current = 10.0 # unspecified units, but does NOT change the potential value. 
    initial_phase = 0.0
#    scaling_frequency = 10221.05558e6 # Angular frequency 
    scaleDef = 2.0 # Arbitrary Ratio of beam width to Synergia length or "Angular Frequency" 
    scaling_frequency = scaleDef*synergia.PH_MKS_c/(2.0*math.pi) # Angular frequency
    part_per_cell = 2 # Result varies little changing this...  as we many cells.. 
    kicks_per_line = 10 # Number of iteration on space charge field kick (irrelevant ? )  
    gridnum = 32 # Tried 32, 48 and 64. Result changes by less < .1 % 
    griddim = (gridnum,gridnum,gridnum)
    num_particles = griddim[0]*griddim[1]*griddim[2] * part_per_cell
    solver = "3D"
    
    print "num_particles =",num_particles
    print "We will use a", solver, "solver"
    
    xwidth=0.0020/scaleDef   # I hope my beam width is 3.2 mm  
    ywidth=xwidth   # same 
    xpwidth=0.0049608  # Irrelevant, we won't tranport the beam 
    ypwidth=0.00768  # same... 
    xOffset = 0.0
    yOffset = 0.0
    rx=0.0 # No offset.   
    dpop = 1.0e-20 # irrelevant as well 
#    bunchLength = xwidth # 3D spherical bunch !... in meters 
    bunchLength = 0.5 # realistic bunch length now..  
#
# Now define the problem in Synergia.  The channel.mad file is from
#  synergia2/examples/channel.  The detail of this beam line are irrelevant, 
# as we will not propagate this bunch in CHEF 
#
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.mad"),"channel",kinetic_energy,
                        scaling_frequency)
			
    gourmet.insert_space_charge_markers(kicks_per_line)

    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=0)
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = -rx)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = rx)
    sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
    # No, set the bunch length to
    sigma_z_meters = xwidth # Make sure we have a 3D isotropic charge. 
    print " sigma_z_meters ", sigma_z_meters
    beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    # Bring it back.. 
    
    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)
    means,sigs = synergia.get_spatial_means_stds(bunch)
    sigXNow = sigs[0]
    sigYNow = sigs[1]
    bunchLength = sigs[2] 
    print "Beam Gamma", beam_parameters.get_gamma()
    print " Bunch sigmas, before propagation, x =  ", \
          sigXNow, " y = ", sigYNow, " z ", bunchLength
    print " Are these in meters, at this stage ???"  
    
    unitsBunch = bunch.get_store().get_units()
    print " Synergia Units X or Y for the bunch ", unitsBunch[0]  
    print " Synergia Units Z for the bunch ", unitsBunch[4]    
    
# How do I get the potential from the bunch after propagation?
# define a new scalar and re-populate with charge from the bunch ?  
# repopulate with charge? This looks like a waste of time !. 
# get the sigma of the bunch.  For now, assume it did not changed by the 
# by the propagation...
# yes, code written by Jim A. and myself.
#      
    sizeNow = (20.0*sigXNow, 20.0*sigXNow, 20.0*sigXNow) # 3D symmetric grid, +- 10 sigma  
# set the offset.. ?? Keep it centered     
    rho = Real_scalar_field(griddim,sizeNow,(0.0,0.0,0.0))
# drop the bunch charge onto the grid..    
    total_charge = deposit_charge_cic(rho,bunch.get_store(),0)
    print " Total Charge ", total_charge # which units? 
    print " Total Current ", bunch.total_current # Which units?
    print " mass  ", bunch.mass # This is in GeV 
# Recompute the potential
#    griddimA = numpy.array([griddim[0], griddim[1], griddim[2]])
    fftwh = Fftw_helper(griddim, False)
    phi = solver_fftw_open(rho,fftwh,0)
#
    iSig=0.0001
    while iSig < 11:
      rr = iSig*sigXNow # at given radius rr .. In what units? 
      loc = numpy.array([rr, rr, rr],'d')
      vv = phi.get_val(loc)/total_charge # Normalize to the total charge (units?) 
#      vvPrime = phi.get_deriv(loc,0)/total_charge # Not available in Python..  
      print " V per Impact Charge ", rr , " = ", vv, " V*r/Q " , vv*rr 
      iSig = iSig+0.5
    
