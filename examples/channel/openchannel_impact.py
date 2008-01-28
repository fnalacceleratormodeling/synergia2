#!/usr/bin/env bwpython

import Numeric
import time
import math


import sys

import synergia
import s2_fish
import UberPkgpy 

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
    newmem0 = synergia.memory()
    print (newmem0-mem0)/math.pow(2,20), "(",newmem0/math.pow(2,20)," total)"
    mem0 = newmem0


if ( __name__ == '__main__'):
    t0 = time.time()
    current = 0.5
    kinetic_energy = 0.0067
    mass = synergia.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6
    part_per_cell = 1
    width_x = 0.004
    pipe_radius = 0.04
    kicks_per_line = 10
    gridnum = int(sys.argv[1])
    griddim = (gridnum,gridnum,gridnum)
    num_particles = griddim[0]*griddim[1]*griddim[2] * part_per_cell
    
    xwidth=0.0012026
    xpwidth=0.0049608
    rx=0.85440
    dpop = 1.0e-20

    ee = synergia.Error_eater()
    ee.start()
    g = synergia.Gourmet("channel.mad","channel",kinetic_energy,
                        scaling_frequency)
    g.insert_space_charge_markers(kicks_per_line)

    bp = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)

    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = xwidth, lam = xpwidth * pz)
    bp.y_params(sigma = xwidth, lam = xpwidth * pz)
    sigma_z_meters = bp.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
    bp.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    bp.correlation_coeffs(xpx = -rx, ypy = rx)

    sys.stdout.flush()
    
    pgrid = synergia.Processor_grid(1)
    cgrid = synergia.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                  "3d open")
    piperad = 0.04
    field = synergia.Field(bp, pgrid, cgrid, piperad)
    b = synergia.Bunch(current, bp, num_particles, pgrid)
    b.generate_particles()
###    b.write_particles_text("oc_particles.h5")
    
    line_length = g.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    diag = synergia.Diagnostics_impact(g.get_initial_u())
    kick_time = 0.0
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
                tk0 = time.time()
                UberPkgpy.Apply_SpaceCharge_external(
                    b.get_beambunch(),
                    pgrid.get_pgrid2d(),
                    field.get_fieldquant(),
                    field.get_compdom(),
                    field.get_period_length(),
                    cgrid.get_bc_num(),
                    field.get_pipe_radius(),
                    tau, 0, scaling_frequency,0)
                tk1 = time.time()
                kick_time += tk1 - tk0
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

    print "elapsed time =",time.time() - t0, "kick time =",kick_time

    #~ diag.write("ocimpact")
    import pylab
    d0 = synergia.Diagnostics_impact_orig("channel0current")

    pylab.plot(d0.s,d0.std[:,synergia.x],'gx',label='Synergia2, no space charge')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    e = synergia.loadfile("envelope_match_channel_0.5A.dat")
    pylab.plot(e[:,0],e[:,1],label='envelope equation')

    pylab.plot(diag.s,diag.std[synergia.x][:],'ro',label='Synergia2 with space charge')

    pylab.show()
