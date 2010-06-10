#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw import solver_fftw_open as solver_fft_open
from s2_solver_fftw import Fftw_helper
from s2_solver_fftw import gather_global_rho
from macro_bunch import get_longitudinal_period_size

from pardebug import pardebug

from mytimer import mytimer
import constraints

import time
import synergia
import sys

from mpi4py import MPI
import Numeric
import MLab
import math

def listify(x):
    if type(x) == type([]) or type(x) == type(()):
        return x
    else:
        return [x]

fftwhs = {}
counter = 0   

def apply_space_charge_kick(shape,size,offset,mbunch_in,tau,
        periodic=True,aperture=None,
                        space_charge=True,
                        impedance=False,
                        pipe_radiusx=None,
                        pipe_radiusy=None,
                        pipe_conduct=None,
                        bunch_spacing=None):
    global counter
    #~ pardebug("start apply_space_charge_kick\n")
    counter += 1
    show_timings=1
    mytimer("misc asck1")
    key = str(shape)
    if not fftwhs.has_key(key):
        fftwhs[key] = Fftw_helper(shape,periodic)
    mbunches = listify(mbunch_in)
    zdensities = []
    xmoms = []
    ymoms = []
    for mbunch in mbunches:
        if aperture:
            constraints.apply_circular_aperture(mbunch.get_store(),aperture)
            mytimer("apply aperture")
        if periodic:
            constraints.apply_longitudinal_periodicity(mbunch.get_store())
            mytimer("apply periodicity")
        mbunch.convert_to_fixedt()
        mytimer("convert")
        if (size == None) or (offset == None):
            means, stds = synergia.get_spatial_means_stds(mbunch)
            if size == None:
                n_sigma = 8.0
                size = list(n_sigma*stds)
            if offset == None:
                offset = list(means)
        if periodic:
            size[2] = get_longitudinal_period_size(mbunch)
            offset[2] = 0
        mytimer("diagnostics")
        if impedance or space_charge:
            #~ print "size =",size
            rho = Real_scalar_field(shape,size,offset)
            total_charge = deposit_charge_cic(rho,mbunch.get_store(),periodic)
            mytimer("deposit")
        if impedance:
            if ((pipe_radiusx == None) or (pipe_radiusy == None) or (pipe_conduct == None)):
                            raise RuntimeError, \
                                "apply_space_charge_kick with impedance != 0.0 requires pipe_radiusx, pipe_radiusy, pipe_conduct and to be specified"
            if ((bunch_spacing == None) and (len(mbunches)>1)):
                            raise RuntimeError, \
                                "apply_space_charge_kick with impedance and multiple bunches requires bunch_spacing to be specified"
            global_rho = Real_scalar_field(shape,size,offset)
            gather_global_rho(rho,global_rho)
            zdensity = Numeric.zeros((shape[2]),'d')
            xmom = Numeric.zeros((shape[2]),'d')
            ymom= Numeric.zeros((shape[2]),'d')
            calculate_rwvars(mbunch.get_store(),zdensity,xmom,ymom,
                     offset[2]-size[2]*0.5,size[2])
            zoffset = 0.0
            rw_kick(rho,zdensity,xmom,ymom,2*tau, 
                    mbunch.get_store(),
                    pipe_radiusx,pipe_radiusy,pipe_conduct,zoffset)
            for index in range(len(zdensities)-1,-1,-1):
                zoffset += bunch_spacing
                rw_kick(rho,zdensities[index],xmoms[index],ymoms[index],2*tau,
                    mbunch.get_store(),
                    pipe_radiusx,pipe_radiusy,pipe_conduct,zoffset)
            zdensities.append(zdensity)
            xmoms.append(xmom)
            ymoms.append(ymom)
                
        if space_charge:
            phi = solver_fft_open(rho,fftwhs[key],periodic,True)
            mytimer("solve")
            full_kick(phi,tau,mbunch.get_store())
            mytimer("full kick")
        mbunch.convert_to_fixedz()
        mytimer("unconvert")
