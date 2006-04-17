#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import physics_constants
import bunch
import diagnostics
import pylab
import beam_parameters
import processor_grid
import computational_grid
import processor_grid
import matching
import time
import field
import math

import sys

#import do_compare_channel
import memory

import UberPkgpy #SpaceChargePkgpy

import apply_map
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

mt = 0
def python_apply_map(map, the_bunch):
    global mt
    t0 = time.time()
    the_bunch.beambunch.Pts1 = \
                             Numeric.matrixmultiply(map,
                                                    the_bunch.beambunch.Pts1)
    mt += time.time() - t0

def wrapped_apply_map(map, the_bunch):
    global mt
    t0 = time.time()
    apply_map.apply_map1(the_bunch.particles(),the_bunch.num_particles_local(),
                         map)
    mt += time.time() - t0

if ( __name__ == '__main__'):
    t0 = time.time()
    current = 0.5
    kinetic_energy = 0.0067
    mass = physics_constants.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6
    part_per_cell = 0.01*100
    width_x = 0.004
    pipe_radius = 0.04
    kicks_per_line = 10
    griddim = (65,17,17)
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,1)
    
    xwidth=0.0012026
    xpwidth=0.0049608
    rx=0.85440
    dpop = 1.0e-20

    printmem()
#    ee = error_eater.Error_eater()
#    ee.start()
    g = gourmet.Gourmet("channel.mad","channel",kinetic_energy,
                        scaling_frequency)
    g.insert_space_charge_markers(2*kicks_per_line)
    g.generate_maps(scaling_frequency)

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)

    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = xwidth, lam = xpwidth * pz)
    bp.y_params(sigma = xwidth, lam = xpwidth * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/scaling_frequency/math.pi
    bp.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    bp.correlation_coeffs(xpx = -rx, ypy = rx)

    print "jfa is here"
    sys.stdout.flush()
    
    printmem("before impact modules")
    pgrid = processor_grid.Processor_grid(1)
    printmem("after pgrid")
    cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                  "trans finite, long periodic round")
    printmem("after cgrid")
    piperad = 0.04
    field = field.Field(bp, pgrid, cgrid, piperad)
    printmem("after field")
    b = bunch.Bunch(current, bp, num_particles, pgrid)
    b.generate_particles()
    printmem("after bunch")

    b.write_particles("initial.dat")
     
    line_length = g.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    b.write_fort(s)

    printmem("before loop")
    for kick in range(0,kicks_per_line):
        if MPI.rank == 0:
            print "-----------------------------------kick %d--------------------------" % kick
        wrapped_apply_map(g.maps[kick*2],b)
        printmem("after leading map %d" % kick)
        sys.stdout.flush()
        UberPkgpy.Apply_SpaceCharge_external(\
            b.get_beambunch(),
            pgrid.get_pgrid2d(),
            field.get_fieldquant(),
            field.get_compdom(),
            field.get_period_length(),
            cgrid.get_bc_num(),
            field.get_pipe_radius(),
            tau, 0, scaling_frequency)
        sys.stdout.flush()
        printmem("after sc kick %d" % kick)
        wrapped_apply_map(g.maps[kick*2+1],b)
        printmem("after trailing map %d" % kick)
        s += line_length/kicks_per_line
        b.write_fort(s)
        printmem("after write fort %d" % kick)
            
    print "elapsed time =",time.time() - t0
    print "map time =",mt

    if MPI.rank == 0:
        import do_compare_channel
        do_compare_channel.doit()
    print "Why does this hang???"
