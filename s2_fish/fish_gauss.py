#!/usr/bin/env python

from GaussSC import *
from mytimer import mytimer

import time
import synergia
import sys

from mpi4py import MPI

def listify(x):
    if type(x) == type([]) or type(x) == type(()):
        return x
    else:
        return [x]

def apply_BasErs_space_charge_kick(mbunch_in,tau,
        space_charge=True,
        impedance=False,
        impedance_pipe_radiusx=None,
        impedance_pipe_radiusy=None,
        pipe_conduct=None,
        bunch_spacing=None):
    show_timings=1
    mytimer("misc asck1")
    mbunches = listify(mbunch_in)
    zdensities = []
    xmoms = []
    ymoms = []
    for mbunch in mbunches:
        if (impedance or space_charge):
            mbunch.convert_to_fixedt()
            mytimer("convert")
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
            n_sigma = 8.0
            means, stds= synergia.get_spatial_means_stds(mbunch)
            ##print "sx, sy, sz =", stds[0], stds[1], stds[2]
            mytimer("diagnostics")
            apply_BasErs_kick(mbunch.get_store(), stds[0], stds[1], tau)
            mytimer("BasErs kick")
        if (impedance or space_charge):        
            mbunch.convert_to_fixedz()
            mytimer("unconvert")
