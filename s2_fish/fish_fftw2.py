#!/usr/bin/env python

from s2_containers import *
from s2_deposit import *
from s2_electric_field import *
from s2_solver_fftw2 import solver_fftw2_open as solver_fft_open

import time
import syn2_diagnostics
import sys

def apply_space_charge_kick(shape,size,offset,mbunch,tau):
	show_timings=0
	t0 = time.time()
	mbunch.convert_to_fixedt()
        t1 = time.time()
	if show_timings:
		print "convert:",t1-t0
 	diag = syn2_diagnostics.Diagnostics([1,1,1,1,1,1])
	coord_stds = diag.get_coord_stds(mbunch)
	n_sigma = 8.0
	size = (coord_stds[0]*n_sigma,coord_stds[1]*n_sigma,coord_stds[2]*n_sigma)
        rho = Real_scalar_field(shape,size,offset)
	t0 = time.time()
        total_charge = deposit_charge_cic(rho,mbunch.get_store(),0)
        t1 = time.time()
	if show_timings:
		print "deposit:",t1-t0
        t0 = time.time()
        phi = solver_fft_open(rho,0)
        t1 = time.time()
	if show_timings:
		print "solve:",t1-t0
        calc_time = 0
        apply_time = 0
	old = 0
	if old:
		for E_axis in range(0,3):
		    t0 = time.time()
		    E = calculate_E_n(phi,E_axis)
		    t1 = time.time()
		    calc_time += t1 -t0
		    t0 = time.time()
		    apply_E_n_kick(E,E_axis,tau,mbunch.get_store())
		    t1 = time.time()
		    apply_time += t1 - t0
		if show_timings:
	 		print "calc E:",calc_time
	 		print "apply:",apply_time
	else:
		t0 = time.time()
		full_kick(phi,tau,mbunch.get_store())
		t1 = time.time()
		full_kick_time = t1-t0;
		if show_timings:
			print "full kick:",full_kick_time
	t0 = time.time()
	mbunch.convert_to_fixedz()
        t1 = time.time()
	if show_timings:
		print "unconvert:",t1-t0

