#!/usr/bin/env python

from s2_solver_cylindrical import *
from macro_bunch import get_longitudinal_period_size
from mytimer import mytimer
import constraints

import time
import synergia
import sys

from mpi4py import MPI
import numpy
import numpy
import math

def listify(x):
    if type(x) == type([]) or type(x) == type(()):
        return x
    else:
        return [x]

counter = 0
def apply_cylindrical_space_charge_kick(shape,radius,mbunch_in,tau,
	periodic=True,
        aperture=None,
        space_charge=True,
        impedance=False,
        impedance_pipe_radiusx=None,
        impedance_pipe_radiusy=None,
        pipe_conduct=None,
        bunch_spacing=None):
    global counter
    counter += 1
    show_timings=1
    mbunches = listify(mbunch_in)
    zdensities = []
    xmoms = []
    ymoms = []
    for mbunch in mbunches:
        if aperture:
            constraints.apply_circular_aperture(mbunch.get_store(),aperture)
            mytimer("apply aperture")
	#if periodic:
        #    constraints.apply_longitudinal_periodicity_t(mbunch.get_store())
        mytimer("apply periodicity")    
        if (impedance or space_charge):
            mbunch.convert_to_fixedt()
	    if periodic:
	        length = get_longitudinal_period_size(mbunch) 
	        constraints.apply_longitudinal_periodicity_z(mbunch.get_store(), length)	    
        if impedance:
            if ((impedance_pipe_radiusx == None) or (impedance_pipe_radiusy == None) or (pipe_conduct == None)):
                            raise RuntimeError, \
                                "apply_space_charge_kick with impedance != 0.0 requires impedance_pipe_radiusx, impedance_pipe_radiusy, pipe_conduct and to be specified"
            if ((bunch_spacing == None) and (len(mbunches)>1)):
                            raise RuntimeError, \
                                "apply_space_charge_kick with impedance and multiple bunches requires bunch_spacing to be specified"
            global_rho = Real_scalar_field(shape,size,offset)
            gather_global_rho(rho,global_rho)
            zdensity = numpy.zeros((shape[2]),'d')
            xmom = numpy.zeros((shape[2]),'d')
            ymom= numpy.zeros((shape[2]),'d')
            calculate_rwvars(mbunch.get_store(),zdensity,xmom,ymom,
                     offset[2]-size[2]*0.5,size[2])
            zoffset = 0.0
            rw_kick(rho,zdensity,xmom,ymom,tau, 
                    mbunch.get_store(),
                    impedance_pipe_radiusx,impedance_pipe_radiusy,pipe_conduct,zoffset)
            for index in range(len(zdensities)-1,-1,-1):
                zoffset += bunch_spacing
                rw_kick(rho,zdensities[index],xmoms[index],ymoms[index],2*tau,
                    mbunch.get_store(),
                    impedance_pipe_radiusx,impedance_pipe_radiusy,pipe_conduct,zoffset)
            zdensities.append(zdensity)
            xmoms.append(xmom)
            ymoms.append(ymom)
        if space_charge:
            coords = numpy.zeros((3,mbunch.local_num),'d')
            get_cylindrical_coords(mbunch.get_store(),coords)
            length = get_longitudinal_period_size(mbunch)
            physical_size = [radius,2*math.pi,length]
            physical_offset = [0.0,0.0,physical_size[2]/2.0]
            periodic = [False,True,True]
            field_domain = Cylindrical_field_domain(radius,length,shape,True)
            rho = numpy.zeros(shape,'d')	   	    
            deposit_charge_cic_cylindrical(field_domain,rho,mbunch.get_store(),coords)		
            phi = numpy.zeros(shape,'d')	   
            solve_cylindrical_finite_periodic(field_domain,rho,phi)  	  
            full_kick_cylindrical(field_domain,phi,tau,mbunch.get_store(),coords)	   
        if (impedance or space_charge):
            mbunch.convert_to_fixedz()
	     
