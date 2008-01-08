#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw import solver_fftw_open as solver_fft_open
from s2_solver_fftw import Fftw_helper
from mytimer import mytimer

import time
import syn2_diagnostics
import sys

from mpi4py import MPI

fftwhs = {}

def apply_space_charge_kick(shape,size,offset,mbunch,tau):
    show_timings=1
    mytimer("misc asck1")
    key = str(shape)
    if not fftwhs.has_key(key):
        fftwhs[key] = Fftw_helper(shape)
    print "jfa: fish_fftw 1"
    print "big one:",mbunch.get_local_particles()[:,0]
    mbunch.convert_to_fixedt()
    print "big two:",mbunch.get_local_particles()[:,0]    
    print "jfa: fish_fftw 2"
    mytimer("convert")
    if (size == None) or (offset == None):
        print "jfa: fish_fftw 3"
        means, stds = syn2_diagnostics.get_spatial_means_stds(mbunch)
        print "jfa: fish_fftw 4"
        print "means,stds =",means,stds
        if size == None:
            n_sigma = 8.0
            size = tuple(n_sigma*stds)
        if offset == None:
            offset = tuple(means)
    print "shape =",shape
    print "size =",size
    print "offset =",offset
    mytimer("diagnostics")
    rho = Real_scalar_field(shape,size,offset)
    print "jfa: fish_fftw 10"
    total_charge = deposit_charge_cic(rho,mbunch.get_store(),0)
    print "jfa: fish_fftw 20"
    mytimer("deposit")
    phi = solver_fft_open(rho,fftwhs[key],0)
    mytimer("solve")
    print "jfa: fish_fftw 30"
    full_kick(phi,tau,mbunch.get_store())
    print "jfa: fish_fftw 40"
    mytimer("full kick")
    mbunch.convert_to_fixedz()
    mytimer("unconvert")
