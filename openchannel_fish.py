#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import physics_constants
import bunch
import diagnostics
import beam_parameters
import processor_grid
import computational_grid
import processor_grid
import matching
import time
import field
import math

import syn2_diagnostics

import fish_fftw as fish
import macro_bunch

import sys

#import do_compare_channel
import memory

import UberPkgpy #SpaceChargePkgpy

import error_eater

from mpi4py import MPI

def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return retval

mem0 = 0.0

def printmem(str=""):
    return None
    global mem0
    if mem0 == 0.0:
        print "memory usage total:",
    else:
        print "memory gain %s:" % str,
    newmem0 = memory.memory()
    print (newmem0-mem0)/math.pow(2,20), "(",newmem0/math.pow(2,20)," total)"
    mem0 = newmem0


if ( __name__ == '__main__'):
    t0 = time.time()
    current = 0.5
    kinetic_energy = 0.0067
    mass = physics_constants.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6
    part_per_cell = 1
    width_x = 0.004
    pipe_radius = 0.04
    kicks_per_line = 10
    gridnum = int(sys.argv[1])
    griddim = (gridnum,gridnum,gridnum)
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,1)
    
    print "num_particles =",num_particles
    xwidth=0.0012026
    xpwidth=0.0049608
    rx=0.85440
    dpop = 1.0e-20

    ee = error_eater.Error_eater()
    ee.start()
    g = gourmet.Gourmet("channel.mad","channel",kinetic_energy,
                        scaling_frequency)
    g.insert_space_charge_markers(kicks_per_line)

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)

    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = xwidth, lam = xpwidth * pz)
    bp.y_params(sigma = xwidth, lam = xpwidth * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/scaling_frequency/math.pi
    bp.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    bp.correlation_coeffs(xpx = -rx, ypy = rx)

    sys.stdout.flush()
    
    pgrid = processor_grid.Processor_grid(1)
    b_impact = bunch.Bunch(current, bp, num_particles, pgrid)
    b_impact.generate_particles()
    
    b = macro_bunch.Macro_bunch(mass,1)
    b.init_from_bunch(b_impact)


    line_length = g.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    diag = syn2_diagnostics.Diagnostics(g.get_initial_u())
    kick_time = 0.0
    for action in g.get_actions():
        if action.is_mapping():
            action.get_data().apply(b.get_local_particles(),
                               b.get_num_particles_local())
            s += action.get_length()
        elif action.is_synergia_action():
            if action.get_synergia_action() == "space charge endpoint":
                diag.add(s,b)
                if not first_action:
                    print "finished space charge kick"
            elif action.get_synergia_action() == "space charge kick":
                tk0 = time.time()
                fish.apply_space_charge_kick(griddim,None,None, b, 2*tau)
                tk1 = time.time()
                kick_time += tk1 - tk0
            elif action.get_synergia_action() == "rfcavity1" or \
                 action.get_synergia_action() == "rfcavity2":
                element = action.get_data()
                u_in = g.get_u(action.get_initial_energy())
                u_out = g.get_u(action.get_final_energy())
                chef_propagate.chef_propagate(
                    b.get_local_particles(), b.get_num_particles_local(),
                    element, action.get_initial_energy(),
                    u_in, u_out)
            else:
                print "unknown action: '%s'" % \
                      action.get_synergia_action()
        else:
            print "action",action.get_type(),"unknown"
        first_action = 0

    print "elapsed time =",time.time() - t0, "kick time =",kick_time

    diag.write("ocfish")
    diag.write_hdf5("ocfish")
    import pylab
    import loadfile

    print \
"""nb: openchannel_fish.py expects to find the results of openchannel_impact.py
in the current directory.
"""
    dimpact = diagnostics.Diagnostics()
    d0 = diagnostics.Diagnostics("channel0current")

    pylab.plot(d0.s,d0.std[:,diagnostics.x],'gx',label='no space charge')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    pylab.plot(dimpact.s,dimpact.std[:,diagnostics.x],'o',label='impact',
               markerfacecolor=None)
    pylab.plot(diag.get_s(),diag.get_stds()[:,syn2_diagnostics.x],'r+',label='barracuda')
    pylab.legend(loc=0)
    pylab.show()
