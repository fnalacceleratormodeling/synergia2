#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import physics_constants
import bunch
import diagnostics
#import pylab
import beam_parameters
import processor_grid
import computational_grid
import processor_grid
import matching
import time
import field
import math

import sys
import memory

import UberPkgpy #SpaceChargePkgpy

import apply_map
import error_eater
import options

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

    myopts = options.Options("syn2booster")
    myopts.add("current",0.035*13,"current",float)
    myopts.add("transverse",1,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("emittance",3.0e-7,"emittance",float)
    myopts.add("zfrac",0.1,"z width as fraction of bucket",float)
    myopts.add("dpop",5.0e-4,"(delta p)/p",float)
    myopts.add("showplot",0,"show plot",int)
    myopts.add("kickspercell",4,"kicks per cell",int)
    myopts.add("plotperiod",4,"update plot every plotperiod steps",int)
    myopts.add("turns",1,"number of booster revolutions",int)
    
    myopts.parse_argv(sys.argv)

### start save
    scaling_frequency = 200.0e6
### end save
    part_per_cell = 4
    pipe_radius = 0.04
    griddim = (33,33,257)
    num_particles = adjust_particles(griddim[0]*griddim[1]*griddim[2] *\
                                     part_per_cell,MPI.size)

    ee = error_eater.Error_eater()
    ee.start()
    g = gourmet.Gourmet("booster_classic.lat","bcel02",0.4,
                        scaling_frequency,myopts.get("maporder"))
    g.insert_space_charge_markers(2*myopts.get("kickspercell"))
    g.generate_mappings()

    bp = beam_parameters.Beam_parameters(physics_constants.PH_NORM_mp,
                                         1.0, 0.4, 0.0, scaling_frequency,
                                         myopts.get("transverse"))

    [sigma_x,sigma_xprime,r_x,\
     sigma_y,sigma_yprime,r_y] = matching.envelope_match(
        myopts.get("emittance"),myopts.get("current"),g)
    
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV
    bp.x_params(sigma = sigma_x, lam = sigma_xprime * pz)
    bp.y_params(sigma = sigma_y, lam = sigma_yprime * pz)
    sigma_z_meters = bp.get_beta()*physics_constants.PH_MKS_c/\
                     scaling_frequency/math.pi * myopts.get("zfrac")
    bp.z_params(sigma = sigma_z_meters, lam = myopts.get("dpop")* pz)
    bp.correlation_coeffs(xpx = -r_x, ypy = -r_y)

    pgrid = processor_grid.Processor_grid(1)
    cgrid = computational_grid.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                  "trans finite, long periodic round")
    piperad = 0.04
    field = field.Field(bp, pgrid, cgrid, piperad)
    b = bunch.Bunch(myopts.get("current"), bp, num_particles, pgrid)
    b.generate_particles()

    b.write_particles("initial.dat")
     
    line_length = g.orbit_length()
    tau = 0.5*line_length/myopts.get("kickspercell")
    s = 0.0
    b.write_fort(s)


    line_x = None
    line_y = None
#     if MPI.rank == 0 and myopts.get("showplot"):
#         pylab.ion()
    steps = 0
    for cell in range(0,myopts.get("turns")*24):
        tcell = time.time()
        for kick in range(0,myopts.get("kickspercell")):
            steps += 1
            if MPI.rank == 0:
                print "cell %d, kick %d" %\
                      (cell,kick)
    ###        wrapped_apply_map(g.maps[kick*2],b)
            g.get_fast_mapping(kick*2).apply(b.particles(), b.num_particles_local())
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
    ###        wrapped_apply_map(g.maps[kick*2+1],b)
            g.get_fast_mapping(kick*2+1).apply(b.particles(), b.num_particles_local())
            s += line_length/myopts.get("kickspercell")
            b.write_fort(s)
            if MPI.rank == 0:
                print "cell time =",time.time() - tcell
#             if MPI.rank == 0 and myopts.get("showplot") and steps % myopts.get("plotperiod") == 0:
#                 d = diagnostics.Diagnostics()
#                 line_x, = pylab.plot(d.s,d.std[:,diagnostics.x],'r-o',label='x')
#                 line_y, = pylab.plot(d.s,d.std[:,diagnostics.y],'b-x',label='y')
#                 pylab.draw()
            
    print "elapsed time =",time.time() - t0

