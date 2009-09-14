#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw import solver_fftw_open as solver_fft_open
from s2_solver_fftw import Fftw_helper
from s2_solver_transverse_fftw import *
from GaussSC import *
from s2_solver_cylindrical import *
#from s2_solver_fftw import gather_global_rho
#from macro_bunch import get_longitudinal_period_size
import constraints

have_impact = False
try:
    import impact
    have_impact = True
except ImportError:
    pass


from mytimer import mytimer
#import constraints

#import time
import synergia
#import sys

from mpi4py import MPI
import numpy
import math


fftwhs = {}

def apply_space_charge_kick_fftw(shape,bunch,tau,periodic,transverse,size=None,offset=None):
	
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
            size[2] = bunch.get_longitudinal_period_size()
            offset[2] = 0
       # mytimer("diagnostics")
	
	#mytimer("deposit")
        rho = Real_scalar_field(shape,size,offset)
        total_charge = deposit_charge_cic(rho,bunch.get_store(),periodic)
        phi = solver_fft_open(rho,fftwhs_key,periodic,True)
        
	#mytimer("solve")
        if transverse:               
           transverse_kick(phi,tau,bunch.get_store(),fftwhs_key,periodic)
        else:
           full_kick(phi,tau,bunch.get_store(),fftwhs_key,periodic)
        #mytimer("full kick")

def apply_space_charge_kick_orbit(shape,mbunch,tau):
       means, stds = synergia.get_spatial_means_stds(mbunch)
       n_sigma = 10.0
       sizes = list(n_sigma*stds)
       offsets = list(means)
       # Recommended minimum grid size 32X32
       nXBins = shape[0]
       nYBins = shape[1]
       eps = 0.001
       includeLocalDensity = True
       solver = TransverseSolver(nXBins,nYBins,eps,includeLocalDensity)
       solver.kick_transverse_charge(mbunch.get_store(), tau,offsets, sizes)
	
def apply_space_charge_kick_gauss(mbunch,tau):
	 n_sigma = 8.0
         means, stds= synergia.get_spatial_means_stds(mbunch)
	 apply_BasErs_kick(mbunch.get_store(), stds[0], stds[1], tau)
	 
def apply_space_charge_kick_cylindrical(shape,radius,mbunch,tau):
    coords = numpy.zeros((3,mbunch.local_num),'d')
    get_cylindrical_coords(mbunch.get_store(),coords)
    length = mbunch.get_longitudinal_period_size()
    physical_size = [radius,2*math.pi,length]
    physical_offset = [0.0,0.0,physical_size[2]/2.0]
    periodic = [False,True,True]
    field_domain = Cylindrical_field_domain(radius,length,shape,True)
    rho = numpy.zeros(shape,'d')	   	    
    deposit_charge_cic_cylindrical(field_domain,rho,mbunch.get_store(),coords)		
    phi = numpy.zeros(shape,'d')	   
    solve_cylindrical_finite_periodic(field_domain,rho,phi)  	  
    full_kick_cylindrical(field_domain,phi,tau,mbunch.get_store(),coords)
	
def apply_space_charge_kick_impact(bunch,tau,pgrid,cgrid,field):		     
	 if ((pgrid == None) or (field == None) or (cgrid == None)):
                        raise RuntimeError, \
                            "propagate with use_impact=True requires pgrid, field and cgrid to be specified"
	 impact.apply_space_charge_kick(
                        bunch.get_beambunch(),
                        pgrid.get_pgrid2d(),
                        field.get_fieldquant(),
                        field.get_compdom(),
                        field.get_period_length(),
                        cgrid.get_bc_num(),
                        field.get_pipe_radius(),
                        tau, 0, bunch.get_scaling_frequency(),0)		    
	 		    
	 		    
			    
# ~ other space charge solver to be put here..(fish_cylindrical, impact...)....
def apply_space_charge_kick(mbunch,space_charge,tau,size=None,offset=None):
    solver=space_charge.get_solver()    		
    if solver=="s2_fish_3d":
        shape=space_charge.get_grid()	    
	   # periodic=space_charge.get_periodic()
        periodic=mbunch.periodic
        transverse=space_charge.get_transverse()
        apply_space_charge_kick_fftw(shape,mbunch,tau,periodic,transverse)
    elif  solver== "s2_fish_gauss": 
        if not(mbunch.periodic):
            raise RuntimeError, "make the bunch periodic and  transverse for gauss solver!"
        apply_space_charge_kick_gauss(mbunch,tau) 
    elif  solver== "s2_fish_cylindrical":
        if not(mbunch.periodic):
            raise RuntimeError, "the cylindrical solver requires a periodic bunch!"
        shape=space_charge.get_grid()
        radius=space_charge.get_radius_cylindrical()
        apply_space_charge_kick_cylindrical(shape,radius,mbunch,tau)	    
    elif solver== "s2_fish_transverse":
        shape = space_charge.get_grid()
        apply_space_charge_kick_orbit(shape,mbunch,tau) 
    elif solver=="impact":
        if not have_impact:
            raise RuntimeError, \
		    "propagate with use_impact=True requires a working impact module"	
        pgrid=space_charge.get_impact_pgrid()
        cgrid=space_charge.get_impact_cgrid()
        field=space_charge.get_impact_field()		    
        apply_space_charge_kick_impact(mbunch,tau,pgrid,cgrid,field)

class SpaceCharge:
    def __init__(self,solver,grid=None, periodic=False,radius_cylindrical=None,transverse=False,impact_pgrid=None, impact_field=None, impact_cgrid=None):
        self.impact_pgrid=impact_pgrid
	self.impact_field=impact_field
	self.impact_cgrid=impact_cgrid
	self.grid=grid		
	self.periodic=periodic
	self.transverse=transverse # only transverse kick
	self.radius_cylindrical=radius_cylindrical
		
	solver_list=["s2_fish_3d", "s2_fish_gauss", "s2_fish_cylindrical", "s2_fish_transverse", "impact"] 
	sv_right=False
	for sv in solver_list:
	    if (solver==sv):
	        sv_right=True
		break	
	if sv_right:
	    self.solver=solver	
	else:
	    raise RuntimeError, \
                            "you chose a wrong solver name or misspeled, see solver_list in space_charge.py "	    
	
    def get_solver(self):
       return self.solver 
    
    
    def get_grid(self):
       return self.grid	
    
    def get_periodic(self):
       return self.periodic
    
    def get_transverse(self):
       return self.transverse
    
    def get_radius_cylindrical(self):
       return self.radius_cylindrical
    
    def get_impact_pgrid(self):
	return self.impact_pgrid 
          
    def get_impact_cgrid(self):
	return self.impact_cgrid
    
    def get_impact_field(self):
	return self.impact_field
     
     
     
