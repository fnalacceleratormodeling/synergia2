#!/usr/bin/env bwpython

import local_paths
import gourmet
import pylab
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

import do_compare
import memory

import SpaceChargePkgpy


def fake_crap(a,b,c,d,e,f,g,h,i,j):
    SpaceChargePkgpy.Apply_SpaceCharge(a,b,c,d,e,f,g,h,i,j)

mem0 = 0.0

def printmem(str=""):
    global mem0
    if mem0 == 0.0:
        print "memory usage total:",
    else:
        print "memory gain %s:" % str,
    newmem0 = memory.memory()
    print (newmem0-mem0)/math.pow(2,20), "(",newmem0/math.pow(2,20)," total)"
    mem0 = newmem0

mt = 0
def apply_map(map, the_bunch):
    global mt
    t0 = time.time()
    the_bunch.beambunch.Pts1 = \
                             Numeric.matrixmultiply(map,
                                                    the_bunch.beambunch.Pts1)
    mt += time.time() - t0

if ( __name__ == '__main__'):
    t0 = time.time()
    current = 2.0
    kinetic_energy = 0.4
    mass = physics_constants.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 200.0e6
    num_particles = 1000000
###    num_particles = 1000
    width_x = 0.004
    pipe_radius = 0.04
    kicks_per_line = 50
    dp_o_p = 0.0
    griddim = (33,33,257)
###    griddim = (17,17,17)
    
    printmem()
    g = gourmet.Gourmet("simplebooster.mad","cell",kinetic_energy)
    g.insert_space_charge_markers(2*kicks_per_line)
    g.generate_maps(scaling_frequency)

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency)
    lf = g.get_lattice_functions()
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV 
    (width_xprime,r_x,emittance) = matching.match_twiss_width(width_x,
                                                              lf.alpha_x[0],
                                                              lf.beta_x[0])

    bp.x_params(sigma = width_x, lam = width_xprime * pz)
    
### use equal horizontal and vertical emittances
    (width_y,width_yprime,r_y) = matching.match_twiss_emittance(emittance,
                                                                lf.alpha_y[0],
                                                                lf.beta_y[0])

    bp.y_params(sigma = width_y, lam = width_yprime * pz)


    bp.z_params(sigma = 0.10, lam = dp_o_p* pz)
    bp.correlation_coeffs(xpx = r_x, ypy = r_y)

    printmem("before impact modules")
    pgrid = processor_grid.Processor_grid()
    printmem("after pgrid")
    cgrid = computational_grid.Computational_grid(33,33,257,
                                                  "trans finite, long periodic round")
    printmem("after cgrid")
    piperad = 0.04
    field = field.Field(bp, pgrid, cgrid, piperad)
    printmem("after field")
    print "<who is printing this?>"
    b = bunch.Bunch(current, bp, num_particles, pgrid)
    print "</who is printing this?>"
    printmem("after bunch")
    
    line_length = g.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    b.write_fort(s)
    printmem("before loop")
    for cell in range(0,2):
        for kick in range(0,kicks_per_line):
            apply_map(g.maps[kick*2],b)
            printmem("after leading map %d" % (kick + kicks_per_line*cell))
            SpaceChargePkgpy.Apply_SpaceCharge_external(b.get_beambunch(),
                                               pgrid.get_pgrid2d(),
                                               field.get_fieldquant(),
                                               field.get_compdom(),
                                               field.get_period_length(),
                                               cgrid.get_bc_num(),
                                               field.get_pipe_radius(),
                                               tau, 0, scaling_frequency)
            printmem("after sc kick %d" % (kick + kicks_per_line*cell))
            apply_map(g.maps[kick*2+1],b)
            printmem("after trailing map %d" % (kick + kicks_per_line*cell))
            s += line_length/kicks_per_line
            b.write_fort(s)
            printmem("after write fort %d" % (kick + kicks_per_line*cell))
            
    print "elapsed time =",time.time() - t0
    print "map time =",mt

    do_compare.doit()
    print "Why does this hang???"
