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
#from s2_deposit import *
from s2_impedance_kick import *


import math
import sys
import numpy

#import string
#import os.path

from mpi4py import MPI

def apply_impedance_kick(bunch,impedance,tau, bunch_index):

    bunchmin, bunchmax = synergia.get_spatial_minmax(bunch)
    rwsize = bunchmax - bunchmin
    impedance_zgrid=impedance.get_z_grid()   
    wake_coeff=numpy.zeros((9), 'd')
    wake_coeff[:]=impedance.get_wake_coeff()[:] # wake_coeff are determined by the pipe_symmetry
    zdensity = numpy.zeros((impedance_zgrid),'d')
    xmom = numpy.zeros((impedance_zgrid),'d')
    ymom= numpy.zeros((impedance_zgrid),'d')
    lnum_part=bunch.local_num
    bin_partition=numpy.zeros((lnum_part),'int')
    calculate_rwvars(bunch.get_store(),zdensity,xmom,ymom,
                     bunchmin[2],rwsize[2], bin_partition)
    
    pipe_radius=impedance.get_pipe_radius()
    pipe_conduct=impedance.get_pipe_conduct()
    wake_factor=impedance.get_wake_factor()
    bool_quad_wake=impedance.get_bool_quad_wake()
    quad_wake_sum=impedance.get_quad_wake_sum()
    cutoff_small_z=impedance.get_cutoff_small_z()
#  quantities needed for beam-beam interaction    
    bunch_spacing=impedance.get_bunch_spacing()
    line_length=impedance.get_orbit_length()
    num_bunches=len(impedance.stored_means)
    nts=len(impedance.stored_means[0].get_means()) # number of turns stored    
    stored_means=numpy.zeros((num_bunches,nts,3),'d')
    stored_buckets=numpy.zeros((num_bunches,nts),'i')
    stored_bunchnp=numpy.zeros((num_bunches,nts),'d')
   
   

    for i in range(0,num_bunches):
        stored_means[i]=impedance.stored_means[i].get_means()        
        stored_buckets[i]=impedance.stored_means[i].get_ibucket()
        stored_bunchnp[i]=impedance.stored_means[i].get_bunchnp()
   
    #the  maximum number of agruments passed to C++ is by defalut 15....grouping parameters
    dparameters=numpy.zeros((7),'d')
    dparameters[0]=rwsize[2]
    dparameters[1]=tau
    dparameters[2]=wake_factor
    dparameters[3]=cutoff_small_z
    dparameters[4]=quad_wake_sum
    dparameters[5]=bunch_spacing
    dparameters[6]=line_length
    
             
    rw_kick(dparameters, bin_partition,  zdensity, xmom, ymom, 
             bunch.get_store(), wake_coeff,  bool_quad_wake, \
             bunch_index, stored_means, stored_buckets,stored_bunchnp)        
   # print "***************************** "

class Impedance:
    '''Defines the parameters and methods necessary to  impedance calculation'''
    def __init__(self, pipe_radius, pipe_conduct, wall_thickness, orbit_length, z_grid, 
        pipe_symmetry="circular",paking_frac=1.0, kick="full",nstored_turns=2):
        self.pipe_radius=pipe_radius
        self.pipe_conduct=pipe_conduct # to agree with eq in chao's book (units s^-1), make pipe_conduct=pipe_conduct/(4*pi*eps_0)	
        self.orbit_length=orbit_length
        self.z_grid=z_grid 
        self.pipe_symmetry=pipe_symmetry
        self.wall_thickness=wall_thickness
        self.paking_frac=paking_frac # <=1.0,  packing fraction is the fraction of elements in the ring  (magnets) which contribute to the res wall impedance 
        self.kick=kick
        self.stored_means=None
        self.nstored_turns=nstored_turns
        self.bunch_spacing=None
        if MPI.COMM_WORLD.Get_rank() == 0:    
                print "impedance kick type=", kick
# pipe symmetry keywords so far, see below  "circular", "x_parallel_plates", "y_parallel_plates", "elliptical"	
        self.r=1.+ wall_thickness*wall_thickness/((wall_thickness+pipe_radius)*(wall_thickness+pipe_radius)) # this definition is for a perfect outer magnet at r=b+t
	#print "geometric factor r=",self.r
#~ cutoff_small is the cut off distance  below which the wake force changes is dependence of z, presenlty we approximate W(z)=0 for z<cutoff_small
        self.cutoff_small_z=4.2*pow(synergia.physics_constants.PH_MKS_c*pipe_radius*pipe_radius*4.*math.pi* \
         synergia.physics_constants.PH_MKS_eps0/(2.*math.pi*pipe_conduct),1./3.)  #~for the factor 4.2 see K. Bane, SLAC-AP-87, 1991
        
 # the depenence W(z) propto 1/sqrt(z) is not valid at distace larger than cutoff_quad where the field penetrates the wall....        
# cutoff_quad is a function of  wall_thickness, still not sure if the form above is the best choice for thick walls (though good for thin walls) 	
	   
	
	
	
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


        if (self.kick == "full") or (self.kick == "transverse"):
            if self.pipe_symmetry=="circular":
    #~ pipe with symmetry regardig pi/2 rotation in xy plane (ex: circular pipe)
                self.bool_quad_wake=False
                ax_dipole=1.0 # please adjust it
                ay_dipole=1.0 # please adjust it
                bx_dipole=0.0 
                by_dipole=0.0 
                a_quad   =0.0
                b_quad   =0.0
                cx_quad  =0.0
                cy_quad  =0.0
            elif self.pipe_symmetry=="x_parallel_plates": # see a. chao prst-ab, 111001, (2002) for parameters
    #~ pipe with (x,z) plane symmetry and (y,z) plane symmetry  (ex: eliptical, parallel plates...)
                self.bool_quad_wake=True
                ax_dipole=math.pi*math.pi/24. #  pi*pi/24 for parallel plates, see a. chao prst-ab, 111001, (2002)
                ay_dipole=math.pi*math.pi/12. # pi*pi/12 for parallel plates, 
                bx_dipole=0.0 
                by_dipole=0.0 
                a_quad  = -math.pi*math.pi/24. #  - pi*pi/24 for parallel plates
                b_quad   =0.0
                cx_quad  =0.0
                cy_quad  =0.0
            elif self.pipe_symmetry=="y_parallel_plates":  # AM guess for the value of the parameters
    #~ pipe with (x,z) plane symmetry and (y,z) plane symmetry  (ex: eliptical, parallel plates...)
                self.bool_quad_wake=True
                ax_dipole=math.pi*math.pi/12. # it's a a guess
                ay_dipole=math.pi*math.pi/24. # it's a a guess
                bx_dipole=0.0 
                by_dipole=0.0 
                a_quad  = math.pi*math.pi/24. #  it's a a guess, not sure about sign
                b_quad   =0.0
                cx_quad  =0.0
                cy_quad  =0.0	    
            elif self.pipe_symmetry=="elliptical": # (consider pipe radius = small axis =b to be along y direction) 
                self.bool_quad_wake=True
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
        elif (self.kick == "longitudinal"):
            self.bool_quad_wake=False
            ax_dipole=0.0 
            ay_dipole=0.0 
            bx_dipole=0.0 
            by_dipole=0.0 
            a_quad   =0.0
            b_quad   =0.0
            cx_quad  =0.0
            cy_quad  =0.0
            
        else:   		
            raise RuntimeError,  "1-- kick can be full, transverse or longitudinal....what did you choose? "
        
        if (self.kick == "full") or (self.kick == "longitudinal") :
            a_monopole=0.25*pipe_radius*pipe_radius
        elif (self.kick == "transverse"):
            a_monopole=0.
        else:          
            raise RuntimeError,  "2-- kick can be full, transverse or longitudinal....what did you choose?"   
                
                
        self.wake_coeff=[ax_dipole,ay_dipole,bx_dipole,by_dipole,a_quad,b_quad,cx_quad,cy_quad, a_monopole]
        
        if (self.bool_quad_wake):
            self.cutoff_quad= int(wall_thickness*wall_thickness*pipe_conduct/ \
                (math.pi*synergia.physics_constants.PH_MKS_eps0*synergia.physics_constants.PH_MKS_c*orbit_length))
            if MPI.COMM_WORLD.Get_rank() == 0:    
                print "long distance cutoff=", self.cutoff_quad ," x  (orbit length =" , self.orbit_length," m )"

            self.quad_wake_sum=self.calculate_quad_wake_sum()
        else: 
            self.cutoff_quad=None
            self.quad_wake_sum=0.
    
        self.wake_factor=paking_frac*synergia.physics_constants.PH_MKS_rp*2.0/ (math.pi*pipe_radius*pipe_radius*pipe_radius)* \
            math.sqrt(4*math.pi*synergia.physics_constants.PH_MKS_eps0*synergia.physics_constants.PH_MKS_c/pipe_conduct)
    
       
     
              
    
    def get_pipe_radius(self):
        return self.pipe_radius
    
    def get_pipe_conduct(self):
        return self.pipe_conduct
    
    def get_wall_thickness(self):
       return self.wall_thickness
    
    def get_orbit_length(self):
        return self.orbit_length
    
    def get_pipe_symmetry(self):
       return self.pipe_symmetry
    
    def get_z_grid(self):
        return self.z_grid
    
    def get_cutoff_small_z(self):
        return self.cutoff_small_z
         
    def get_wake_coeff(self):
        return self.wake_coeff
    
    def get_bool_quad_wake(self):
        return self.bool_quad_wake
    
    def get_quad_wake_sum(self):
         return self.quad_wake_sum
    
    def get_wake_factor(self):
         return self.wake_factor
    
    
    def calculate_quad_wake_sum(self):
        orbit_length= self.orbit_length
        cutoff=self.cutoff_quad	
        t=self.wall_thickness
        sigma=self.pipe_conduct
        quad_wake_sum=0.
        for turn in range(self.nstored_turns,cutoff+1):		 
            quad_wake_sum +=1./math.sqrt(turn) 
             	    	    
        quad_wake_sum *= 1./math.sqrt(orbit_length)
	#print " quad_wake_sum  until cut off=", quad_wake_sum
# quad_wake_exp is the exponentially decaying term at  distance > cutoff*orbit_length
        r=self.r
        b=self.pipe_radius 
#  W(z) is prop to exp(-aa*z) at large z
        aa=(synergia.physics_constants.PH_MKS_c*synergia.physics_constants.PH_MKS_eps0)/(sigma*t*r*b) 
        quad_wake_exp=math.exp(-aa*(cutoff+1)*orbit_length)/(orbit_length*aa)
        quad_wake_exp *= math.sqrt(math.pi*synergia.physics_constants.PH_MKS_c*synergia.physics_constants.PH_MKS_eps0/sigma)/t 

	#print "quad_wake_exp (after cutoff contribution) =", quad_wake_exp
        quad_wake_sum += quad_wake_exp
	#print " quad_wake_sum =", quad_wake_sum
        return  quad_wake_sum   
	

    
    def get_means(self):
        return numpy.array(self.means)
    
    def get_bunch_spacing(self):
        return self.bunch_spacing 
        