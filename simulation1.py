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
import matching
import time

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
    current = 0
    kinetic_energy = 0.4
    mass = physics_constants.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 200.0e6
    num_particles = 1000
    width_x = 0.004

    kicks_per_line = 4
    g = gourmet.Gourmet("simplebooster.mad","cell",kinetic_energy)
    g.insert_space_charge_markers(kicks_per_line)
    g.generate_maps(200.0e6)

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency)
    lf = g.get_lattice_functions()
    pz = bp.gamma() * bp.beta() * bp.mass_GeV 
    (width_xprime,r_x,emittance) = matching.match_twiss_width(width_x,
                                                              lf.alpha_x[0],
                                                              lf.beta_x[0])

    bp.x_params(sigma = width_x, lam = width_xprime * pz)
    
### use equal horizontal and vertical emittances
    (width_y,width_yprime,r_y) = matching.match_twiss_emittance(emittance,
                                                                lf.alpha_y[0],
                                                                lf.beta_y[0])

    bp.y_params(sigma = width_y, lam = width_yprime * pz)

    dp_o_p = 0.0 # myopts.get_value("dpop")
    bp.z_params(sigma = 0.10, lam = dp_o_p* pz)
    bp.correlation_coeffs(xpx = r_x, ypy = r_y)

    print "<who is printing this?>"
    b = bunch.Bunch(current, bp, num_particles)
    print "</who is printing this?>"
    line_length = g.orbit_length()
    s = 0.0
    b.write_fort(s)
    for cell in range(0,24):
        for map in g.maps:
            apply_map(map, b)
            s += line_length/kicks_per_line
            b.write_fort(s)
            
    print "elapsed time =",time.time() - t0
    print "map time =",mt
    d = diagnostics.Diagnostics()
    dold = diagnostics.Diagnostics("/home3/amundson/work/Layer-head/transition/cl.01")

    pylab.plot(d.s,d.std[:,diagnostics.x])
    pylab.plot(dold.s,dold.std[:,diagnostics.x])
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')
    pylab.figure()
    pylab.plot(d.s,d.corr[:,diagnostics.x_xprime])    
    pylab.plot(dold.s,dold.corr[:,diagnostics.x_xprime])    
    pylab.ylabel("r_xx' (unitless)")
    pylab.xlabel('s (m)')
    pylab.show()
    print "Why does this hang???"
