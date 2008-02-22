#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw import solver_fftw_open as solver_fft_open
from s2_solver_fftw import Fftw_helper
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

def get_longitudinal_period_size(mbunch):
    mbs = mbunch.get_store()
    ref_particle = mbs.get_ref_particle()
    units = mbs.get_units()
    gamma = -ref_particle[5]
    beta = math.sqrt(1.0-1.0/gamma**2)
    size = 2.0*math.pi*gamma*beta/units[0]
    return size

def apply_space_charge_kick(shape,size,offset,mbunch,tau,
        periodic=True,aperture=None):
    global counter
    counter += 1
    show_timings=1
    mytimer("misc asck1")
    key = str(shape)
    if not fftwhs.has_key(key):
        fftwhs[key] = Fftw_helper(shape)
    if aperture:
        constraints.apply_circular_aperture(mbunch.get_store(),aperture)
        mytimer("apply aperture")
    if periodic:
        constraints.apply_longitudinal_periodicity(mbunch.get_store())
        mytimer("apply periodicity")
    mbunch.convert_to_fixedt()
    mbunch.write_particles("bunch_fixedt_%02d.dat" % counter)
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
    rho.write_to_file("rho_%02d.dat"%counter)
    mytimer("deposit")
    phi = solver_fft_open(rho,fftwhs[key],periodic)
    phi.write_to_file("phi_%02d.dat"%counter)
    mytimer("solve")
    full_kick(phi,tau,mbunch.get_store())
    mytimer("full kick")
    mbunch.convert_to_fixedz()
    mytimer("unconvert")
