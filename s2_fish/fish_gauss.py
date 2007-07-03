#!/usr/bin/env python

from GaussSC import *
from mytimer import mytimer

import time
import syn2_diagnostics
import sys

from mpi4py import MPI


def apply_BasErs_space_charge_kick(mbunch,tau, PartPersigmaZ):
    show_timings=1
    mytimer("misc asck1")

    mbunch.convert_to_fixedt()
    mytimer("convert")

    n_sigma = 8.0
    means, stds= syn2_diagnostics.get_spatial_means_stds(mbunch)
    mytimer("diagnostics")

    apply_BasErs_kick(mbunch.get_store(), stds[0], stds[1], stds[2],tau, PartPersigmaZ)
    mytimer("BasErs kick")

    mbunch.convert_to_fixedz()
    mytimer("unconvert")
