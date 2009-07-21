#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw import solver_fftw_open as solver_fft_open
from s2_solver_fftw import Fftw_helper
#from s2_solver_fftw import gather_global_rho
from macro_bunch import get_longitudinal_period_size




from mytimer import mytimer
#import constraints

#import time
import synergia
#import sys

from mpi4py import MPI
#import numpy
#import math


fftwhs = {}

def apply_space_charge_kick_fftw(shape,size,offset,bunch,tau,periodic,transverse):
	
    key = str(shape)
    if not fftwhs.has_key(key):	    
        fftwhs[key] = Fftw_helper(shape,periodic)
    fftwhs_key=fftwhs[key]
	
    if (size == None) or (offset == None):
        means, stds = synergia.get_spatial_means_stds(bunch)
        
	if size == None:
             n_sigma = 10.0
             size = list(n_sigma*stds)

        if offset == None:
            offset = list(means)
 
        if periodic:
            size[2] = get_longitudinal_period_size(bunch)
            offset[2] = 0
       # mytimer("diagnostics")
	
	#mytimer("deposit")
        rho = Real_scalar_field(shape,size,offset)
        total_charge = deposit_charge_cic(rho,bunch.get_store(),periodic)
        phi = solver_fft_open(rho,fftwhs_key,periodic,True)
        
	#mytimer("solve")
        if transverse:               
           transverse_kick(phi,tau,mbunch.get_store(),fftwhs_key,periodic)
        else:
           full_kick(phi,tau,bunch.get_store(),fftwhs_key,periodic)
        #mytimer("full kick")


# ~ other space charge solver to be put here..(fish_cylindrical, fish_gauss etc...)....