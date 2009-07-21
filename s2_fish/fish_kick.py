#!/usr/bin/env python

#from s2_containers import *
#from s2_deposit import *
#from s2_electric_field import *
#from s2_solver_fftw import solver_fftw_open as solver_fft_open
#from s2_solver_fftw import Fftw_helper
#from s2_solver_fftw import gather_global_rho
from macro_bunch import get_longitudinal_period_size
from impedance import *
from space_charge import *

from pardebug import pardebug

from mytimer import mytimer
import constraints

import time
import synergia
import sys

from mpi4py import MPI
#import numpy
#import math





def listify(x):
    if type(x) == type([]) or type(x) == type(()):
        return x
    else:
        return [x]


counter = 0   

def apply_kick(shape,size,offset,mbunch_in,tau,
        periodic=True,aperture=None,
                        space_charge=True,
                        impedance=None,                        
                        transverse=False):

    
    global counter
    counter += 1
    show_timings=1
         
    mbunches = listify(mbunch_in) 
    for mbunch in mbunches:
        if aperture:
            constraints.apply_circular_aperture(mbunch.get_store(),aperture)
            mytimer("apply aperture")  
        mbunch.convert_to_fixedt()
	mytimer("convert")
	
	if periodic:	
	    length=get_longitudinal_period_size(mbunch)	
	    constraints.apply_longitudinal_periodicity_z(mbunch.get_store(),length)
	    mytimer("apply periodicity")
        	
		

	if space_charge:
	    apply_space_charge_kick_fftw(shape,size,offset,mbunch,tau,periodic,transverse)	
           
  
        if impedance:
	    if (len(mbunches) >1): 
                print "impedance for multiple bunches not implemented yet"
	        sys.exit(1)
	    apply_impedance_kick(mbunch,impedance,tau) 
	
	mbunch.convert_to_fixedz()
        mytimer("unconvert")
