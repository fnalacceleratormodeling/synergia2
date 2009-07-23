#!/usr/bin/env python


#from bmlfactory import *
#from basic_toolkit import *
#from mxyzptlk import *
#from beamline import *
#from physics_toolkit import *
#from physics_constants import *
#import mappers
#import chef_propagate

import synergia
from s2_deposit import *
from s2_electric_field import *


import math
#import sys
import numpy
#import string
#import os.path

#from mpi4py import MPI

def apply_impedance_kick(bunch,impedance,tau):
    bunchmin, bunchmax = synergia.get_spatial_minmax(bunch)
    rwsize = bunchmax - bunchmin
    impedance_zgrid=impedance.get_z_grid()   
    wake_coeff=numpy.zeros((8), 'd')
    wake_coeff[:]=impedance.get_wake_coeff()[:]
    zdensity = numpy.zeros((impedance_zgrid),'d')
    xmom = numpy.zeros((impedance_zgrid),'d')
    ymom= numpy.zeros((impedance_zgrid),'d')
    lnum_part=bunch.local_num
    bin_partition=numpy.zeros((lnum_part),'int')
    calculate_rwvars(bunch.get_store(),zdensity,xmom,ymom,
                     bunchmin[2],rwsize[2], bin_partition)
    
    
    pipe_radius=impedance.get_pipe_radius()
    pipe_conduct=impedance.get_pipe_conduct()
    quad_wake=impedance.get_quad_wake()
    orbit_length=impedance.get_orbit_length()
    quad_wake_sum=impedance.get_quad_wake_sum()
    cutoff_small_z=impedance.get_cutoff_small_z()
    rw_kick(bunchmin[2],rwsize[2],bin_partition,
                    zdensity,xmom,ymom,tau, 
                    bunch.get_store(),
                    pipe_radius,pipe_conduct,cutoff_small_z,wake_coeff,orbit_length, quad_wake_sum, quad_wake)
    
    

class Impedance:
    '''Defines the parameters and methods necessary to  impedance calculation'''
    def __init__(self, pipe_radius, pipe_conduct, orbit_length, z_grid, pipe_symmetry="circular"):
	    
	self.pipe_radius=pipe_radius
	self.pipe_conduct=pipe_conduct # to agree with eq in chao's book (units s^-1), make pipe_conduct=pipe_conduct/(4*pi*eps_0)
	self.orbit_length=orbit_length
	self.z_grid=z_grid 
	self.cutoff_small_z=4.2*pow(synergia.physics_constants.PH_MKS_c*pipe_radius*pipe_radius*4.*math.pi* synergia.physics_constants.PH_MKS_eps0/(2.*math.pi*pipe_conduct),1./3.)
	print "cutoff_small_z=",self.cutoff_small_z #~for the factor 4.2 see K. Bane, SLAC-AP-87, 1991
	self.pipe_symmetry=pipe_symmetry
# pipe symmetry keywords so far, see below  "circular", "x parallel plates", "y parallel plates", "elliptical"	
        self.wall_thickness=None
	self.quad_wake=None
        self.wake_coeff=None 
	self.cutoff_quad=1000 # cutoff shoul be calculated as a function of  wall_thickness
	
	
	
	
	
	
#see S. Heifets, SLAC/AP110, January 1998....for the following 

       #xkick=wake_factor*(ax_dipole*x_leading+
                     #bx_dipole*y_leading+
                    #(cx_quad+a_quad*x_test+b_quad*y_test));
 
       #ykick=wake_factor*(ay_dipole*y_leading+
                    #by_dipole*x_leading+
                    #(cy_quad-a_quad*y_test+b_quad*x_test));

#some parameters are zero due to symmetry

#general input parameter (ax_dipole, ay_dipole, bx_dipole, by_dipole,a_quad, b_quad, cx_quad, cy_quad)

    #~if the pipe has (x,z) plane symmetry 
    #~bx_dipole=by_dipole=0
    #~b_quad=0
    #~cy_quad =0


   #~if the pipe has (y,z) plane symmetry 
   #~bx_dipole=by_dipole=0
   #~b_quad=0
   #~cx_quad =0



        if self.pipe_symmetry=="circular":
#~ pipe with symmetry regardig pi/2 rotation in xy plane (ex: circular pipe)
            self.quad_wake=False
	    ax_dipole=1.0 # please adjust it
	    ay_dipole=1.0 # please adjust it
	    bx_dipole=0.0 
            by_dipole=0.0 
            a_quad   =0.0
            b_quad   =0.0
            cx_quad  =0.0
            cy_quad  =0.0
	elif self.pipe_symmetry=="x parallel plates": # see a. chao prst-ab, 111001, (2002) for parameters
#~ pipe with (x,z) plane symmetry and (y,z) plane symmetry  (ex: eliptical, parallel plates...)
	    self.quad_wake=True
            ax_dipole=math.pi*math.pi/24. #  pi*pi/24 for parallel plates, see a. chao prst-ab, 111001, (2002)
	    ay_dipole=math.pi*math.pi/12. # pi*pi/12 for parallel plates, 
	    bx_dipole=0.0 
            by_dipole=0.0 
            a_quad  = -math.pi*math.pi/24. #  - pi*pi/24 for parallel plates
            b_quad   =0.0
            cx_quad  =0.0
            cy_quad  =0.0
        elif self.pipe_symmetry=="y parallel plates":  # AM guess for the value of the parameters
#~ pipe with (x,z) plane symmetry and (y,z) plane symmetry  (ex: eliptical, parallel plates...)
	    self.quad_wake=True
            ax_dipole=math.pi*math.pi/12. # it's a a guess
	    ay_dipole=math.pi*math.pi/24. # it's a a guess
	    bx_dipole=0.0 
            by_dipole=0.0 
            a_quad  = math.pi*math.pi/24. #  it's a a guess, not sure about sign
            b_quad   =0.0
            cx_quad  =0.0
            cy_quad  =0.0	    
	elif self.pipe_symmetry=="elliptical": # (consider pipe radius = small axis =b to be along y direction) 
	    self.quad_wake=True
            ax_dipole=math.pi*math.pi/24. # please adjust it, between (pi*pi/24,1) ?? am guess
            ay_dipole=math.pi*math.pi/12. # please adjust it, between (pi*pi/12,1) ?? am guess
	    bx_dipole=0.0 
            by_dipole=0.0 
            a_quad  = -math.pi*math.pi/24. # please adjust it, between (- pi*pi/24, 0) ?? am guess
            b_quad   =0.0
            cx_quad  =0.0
            cy_quad  =0.0
	else:
	    raise RuntimeError,  "pipe symmetry wrongly specified"			
	self.wake_coeff=[ax_dipole,ay_dipole,bx_dipole,by_dipole,a_quad,b_quad,cx_quad,cy_quad]
        
	if (self.quad_wake):
	    self.quad_wake_sum=self.calculate_quad_wake_sum()
	else: 
            self.quad_wake_sum=0.
    
    
    def get_pipe_radius(self):
        return self.pipe_radius
    
    def get_pipe_conduct(self):
        return self.pipe_conduct
    
    def get_orbit_length(self):
        return self.orbit_length
    
    def get_z_grid(self):
        return self.z_grid
    
    def get_cutoff_small_z(self):
	 return self.cutoff_small_z
         
    def get_wake_coeff(self):
	   return self.wake_coeff
    
    def get_quad_wake(self):
	   return self.quad_wake
    
    def get_quad_wake_sum(self):
         return self.quad_wake_sum
    
    def calculate_quad_wake_sum(self):
	orbit_length= self.orbit_length
	cutoff=self.cutoff_quad
	quad_wake_sum=0.
	for turn in range(1,cutoff+1):		 
	    quad_wake_sum +=1./math.sqrt(turn*orbit_length)  	    
	
	quad_wake_exp=0.0	# this is the exponentially decaying term at  distance > cutoff*orbit_length    
	
	quad_wake_sum += quad_wake_exp
	return  quad_wake_sum   
	
   	