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

if ( __name__ == '__main__'):
    t0 = time.time()
    current = 1.5
    kinetic_energy = 0.5
    mass = physics_constants.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 1.0e8
    part_per_cell = 1
    width_x = 0.004
    kicks_per_line = 50
    griddim = (65,17,17)
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,1)
    

    ee = error_eater.Error_eater()
    ee.start()
    g = gourmet.Gourmet("fodo.lat","fodo",kinetic_energy,
                        scaling_frequency)

    if current < 1.0e-6:
        print "using zero-current matching"
        (alpha_x, alpha_y, beta_x, beta_y) = matching.get_alpha_beta(g)
        print "(alpha_x, alpha_y, beta_x, beta_y) =",(alpha_x, alpha_y, beta_x, beta_y)
        (xpwidth,rx,emittance) = matching.match_twiss_width(width_x,alpha_x,beta_x)
        xwidth=width_x
        (ywidth,ypwidth,ry)= matching.match_twiss_emittance(emittance,
                alpha_y,beta_y)
        print "(xwidth,xpwidth,rx,ywidth,ypwidth,ry) =",(xwidth,xpwidth,rx,ywidth,ypwidth,ry)
    else:
        print "using envelope matching"
        # get emittance as though we were doing zero-current matching
        (alpha_x, alpha_y, beta_x, beta_y) = matching.get_alpha_beta(g)
        (xpwidth,rx,emittance) = matching.match_twiss_width(width_x,alpha_x,beta_x)
        (xwidth,xpwidth,rx,ywidth,ypwidth,ry) = \
            matching.envelope_match(emittance,current,g)
        rx = -rx
        ry = -ry
        print "(xwidth,xpwidth,rx,ywidth,ypwidth,ry) =",(xwidth,xpwidth,rx,ywidth,ypwidth,ry)

    dpop = 1.0e-4

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)

    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = xwidth, lam = xpwidth * pz)
    bp.y_params(sigma = ywidth, lam = ypwidth * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/scaling_frequency/math.pi
    bp.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    bp.correlation_coeffs(xpx = rx, ypy = ry)

    g.insert_space_charge_markers(kicks_per_line)

    sys.stdout.flush()
    
    pgrid = processor_grid.Processor_grid(1)
    cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                  "trans finite, long periodic round")
    piperad = 0.04
    field = field.Field(bp, pgrid, cgrid, piperad)
    b = bunch.Bunch(current, bp, num_particles, pgrid)
    b.generate_particles()

    line_length = g.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    for cell in range(0,2):
        for action in g.get_actions():
            if action.is_mapping():
                action.get_data().apply(b.particles(),
                                   b.num_particles_local())
                s += action.get_length()
            elif action.is_synergia_action():
                if action.get_synergia_action() == "space charge endpoint":
                    b.write_fort(s)
                    #~ if not first_action:
                        #~ print "finished space charge kick"
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

    if MPI.COMM_WORLD.Get_rank() == 0:
        import pylab
        d = diagnostics.Diagnostics()
        pylab.plot(d.s,d.std[:,diagnostics.x],'b',label='x')
        pylab.plot(d.s,d.std[:,diagnostics.y],'r',label='y')
        pylab.xlabel('s (m)')
        #~ pylab.ylabel('width (m)')
        #~ d = diagnostics.Diagnostics("wtf_pipe08")
        #~ pylab.plot(d.s,d.std[:,diagnostics.x],'g',label='x')
        #~ pylab.plot(d.s,d.std[:,diagnostics.y],'y',label='y')
        
        pylab.show()
