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
    if not fftwhs.has_key(shape):
        fftwhs[shape] = Fftw_helper(shape)
    mbunch.convert_to_fixedt()
    mytimer("convert")
    n_sigma = 8.0
    stdx,stdy,stdz = syn2_diagnostics.get_spatial_stds(mbunch)
    size = (stdx*n_sigma,stdy*n_sigma,stdz*n_sigma)
    mytimer("diagnostics")
    rho = Real_scalar_field(shape,size,offset)
    total_charge = deposit_charge_cic(rho,mbunch.get_store(),0)
    mytimer("deposit")
    phi = solver_fft_open(rho,fftwhs[shape],0)
    mytimer("solve")
    old = 0
    if old:
        for E_axis in range(0,3):
            E = calculate_E_n(phi,E_axis)
            calc_time += t1 -t0
            apply_E_n_kick(E,E_axis,tau,mbunch.get_store())
        mytimer("full kick")
    else:
        full_kick(phi,tau,mbunch.get_store())
        mytimer("full kick")
    mbunch.convert_to_fixedz()
    mytimer("unconvert")
