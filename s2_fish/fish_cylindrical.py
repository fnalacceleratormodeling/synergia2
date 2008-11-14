#!/usr/bin/env python

from s2_solver_cylindrical import *
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

counter = 0

def apply_cylindrical_space_charge_kick(shape,radius,mbunch,tau,
        aperture=None):
    global counter
    counter += 1
    show_timings=1
    if aperture:
        constraints.apply_circular_aperture(mbunch.get_store(),aperture)
        mytimer("apply aperture")
    mbunch.convert_to_fixedt()
    coords = Numeric.zeros((3,mbunch.local_num),'d')
    get_cylindrical_coords(mbunch.get_store(),coords)
    length = get_longitudinal_period_size(mbunch)
    physical_size = [radius,2*math.pi,length]
    physical_offset = [0.0,0.0,physical_size[2]/2.0]
    periodic = [False,True,True]
    field_domain = Field_domain(physical_size,physical_offset,shape,periodic)
    rho = Numeric.zeros(shape,'d')
    print "jfa: about to deposit_charge_cic_cylindrical"
    deposit_charge_cic_cylindrical(field_domain,rho,mbunch.get_store(),coords)
    print "rho =",rho
    phi = Numeric.zeros(shape,'d')
    solve_cylindrical_finite_periodic(field_domain,rho,phi)

    full_kick_cylindrical(field_domain,phi,tau,mbunch.get_store(),coords)
    mbunch.convert_to_fixedz()
