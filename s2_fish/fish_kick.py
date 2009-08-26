#!/usr/bin/env python

from macro_bunch import get_longitudinal_period_size
from impedance import *
from space_charge import *
from s2_containers import *
from s2_solver_transverse_fftw import *
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

def apply_kick(mbunch_in,tau, space_charge=None,impedance=None, periodic_bunch=False, aperture=None):
    
    
    global counter 
    counter += 1
    show_timings=1
    # ONLY FOR DEBUG PURPOSES
    #if counter > 2:
    #	return
       
    mbunches = listify(mbunch_in) 
    for mbunch in mbunches:    
        if aperture:
            constraints.apply_circular_aperture(mbunch.get_store(),aperture)
            mytimer("apply aperture")  
        mbunch.convert_to_fixedt()
	mytimer("convert")
	
	if periodic_bunch:	
	    length=get_longitudinal_period_size(mbunch)	
	    constraints.apply_longitudinal_periodicity_z(mbunch.get_store(),length)
	    mytimer("apply periodicity")
        
	    
	if space_charge:
	    apply_space_charge_kick(mbunch,space_charge,tau)   
	    
        if impedance:
	    if (len(mbunches) >1): 
                print "impedance for multiple bunches not implemented yet"
	        sys.exit(1)
	    apply_impedance_kick(mbunch,impedance,tau) 
	
	mbunch.convert_to_fixedz()
        mytimer("unconvert")


	
