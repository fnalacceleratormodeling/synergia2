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
    show_timings=0
    mytimer("misc asck1")
    key = str(shape)
    if not fftwhs.has_key(key):
        fftwhs[key] = Fftw_helper(shape)
    mbunch.convert_to_fixedt()
    mytimer("convert")
    if (size == None) or (offset == None):
        means, stds = syn2_diagnostics.get_spatial_means_stds(mbunch)
        if size == None:
            n_sigma = 8.0
            size = tuple(n_sigma*stds)
        if offset == None:
            offset = tuple(means)
    mytimer("diagnostics")
    rho = Real_scalar_field(shape,size,offset)
    total_charge = deposit_charge_cic(rho,mbunch.get_store(),0)
    mytimer("deposit")
    phi = solver_fft_open(rho,fftwhs[key],0)
    mytimer("solve")
    full_kick(phi,tau,mbunch.get_store())
    mytimer("full kick")
    mbunch.convert_to_fixedz()
    mytimer("unconvert")
