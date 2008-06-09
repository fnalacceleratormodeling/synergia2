#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw import solver_fftw_open as solver_fft_open
from s2_solver_fftw import Fftw_helper
from macro_bunch import get_longitudinal_period_size

from mytimer import mytimer
import constraints

import time
import synergia
import sys

from mpi4py import MPI
import Numeric
import MLab
import math

fftwhs = {}
counter = 0   

def apply_space_charge_kick(shape,size,offset,mbunch,tau,
        periodic=True,aperture=None):
    global counter
    counter += 1
    show_timings=1
    mytimer("misc asck1")
    key = str(shape)
    if not fftwhs.has_key(key):
        fftwhs[key] = Fftw_helper(shape,periodic)
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
    rho = Real_scalar_field(shape,size,offset)
    total_charge = deposit_charge_cic(rho,mbunch.get_store(),periodic)
    rho.get_points().print_("rho")
    mytimer("deposit")
    phi = solver_fft_open(rho,fftwhs[key],periodic)
    mytimer("solve")
    full_kick(phi,tau,mbunch.get_store())
    mytimer("full kick")
    mbunch.convert_to_fixedz()
    mytimer("unconvert")
