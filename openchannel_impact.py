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

import syn2_diagnostics_impact

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
    part_per_cell = 0.01*100
    width_x = 0.004
    pipe_radius = 0.04
    kicks_per_line = 10
    griddim = (64,64,64)
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,1)
    
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
    cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                  "3d open")
    piperad = 0.04
    field = field.Field(bp, pgrid, cgrid, piperad)
    b = bunch.Bunch(current, bp, num_particles, pgrid)
    b.generate_particles()
###    b.write_particles_text("oc_particles.h5")
    
    line_length = g.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    diag = syn2_diagnostics_impact.Diagnostics(g.get_initial_u())
    for action in g.get_actions():
        if action.is_mapping():
            action.get_data().apply(b.particles(),
                               b.num_particles_local())
            s += action.get_length()
        elif action.is_synergia_action():
            if action.get_synergia_action() == "space charge endpoint":
                diag.add(s,b)
                b.write_fort(s)
                if not first_action:
                    print "finished space charge kick"
            elif action.get_synergia_action() == "space charge kick":
                UberPkgpy.Apply_SpaceCharge_external(
                    b.get_beambunch(),
                    pgrid.get_pgrid2d(),
                    field.get_fieldquant(),
                    field.get_compdom(),
                    field.get_period_length(),
                    cgrid.get_bc_num(),
                    field.get_pipe_radius(),
                    tau, 0, scaling_frequency,0)
            elif action.get_synergia_action() == "rfcavity1" or \
                 action.get_synergia_action() == "rfcavity2":
                element = action.get_data()
                u_in = g.get_u(action.get_initial_energy())
                u_out = g.get_u(action.get_final_energy())
                chef_propagate.chef_propagate(
                    b.particles(), b.num_particles_local(),
                    element, action.get_initial_energy(),
                    u_in, u_out)
            else:
                print "unknown action: '%s'" % \
                      action.get_synergia_action()
        else:
            print "action",action.get_type(),"unknown"
        first_action = 0

    print "elapsed time =",time.time() - t0

    diag.write("ocimpact")
    import pylab
    import loadfile

    d = diagnostics.Diagnostics()

    d0 = diagnostics.Diagnostics("channel0current")

    pylab.plot(d0.s,d0.std[:,diagnostics.x],'gx',label='Synergia2, no space charge')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    e = loadfile.loadfile("envelope_match_channel_0.5A.dat")
    pylab.plot(e[:,0],e[:,1],label='envelope equation')

    dold = diagnostics.Diagnostics(".")
#    pylab.plot(dold.s,dold.std[:,diagnostics.x],'yo')

    pylab.plot(d.s,d.std[:,diagnostics.x],'ro',label='Synergia2 with space charge')
    pylab.plot(diag.s,diag.std[syn2_diagnostics_impact.x][:],'b+')
###    pylab.legend(loc=0)

    pylab.show()
