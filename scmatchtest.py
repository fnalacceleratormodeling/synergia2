#!/usr/bin/env bwpython

import local_paths
import gourmet
import Numeric
import LinearAlgebra
import MLab
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
from math import cos, sin

import sys

import pylab

#import plot_bunch
import memory

import UberPkgpy #SpaceChargePkgpy

import mappers
import error_eater

import octapy
import function_cache

from mpi4py import MPI

def adjust_particles(base,procs):
    retval = base
    multiple = base/(procs * 10.0)
    if not multiple == round(multiple):
        retval = round(multiple + 1) * \
                   (procs * 10)
    return retval

def get_alpha_beta(my_gourmet):
    mymap = my_gourmet.get_single_linear_map()
    u = my_gourmet.get_u()
    mxx = mymap[0,0]
    mxpxp = mymap[1,1]
    mxxp = mymap[0,1]
    cos_mu = (mxx+mxpxp)/2.0
    mu = math.acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if mxxp/sin(mu) < 0:
        mu = 2*math.pi - mu
    beta_x = mxxp/sin(mu)*u[1]/u[0]
    alpha_x = (mxx-mxpxp)/(2.0*sin(mu))

    myy = mymap[2,2]
    mypyp = mymap[3,3]
    myyp = mymap[2,3]
    cos_mu = (myy+mypyp)/2.0
    mu = math.acos(cos_mu)
    # beta function is positive
    # use this to pick branch
    if myyp/sin(mu) < 0:
        mu = 2*math.pi - mu
    beta_y = myyp/sin(mu)*u[3]/u[2]
    alpha_y = (myy-mypyp)/(2.0*sin(mu))

    return (alpha_x, alpha_y, beta_x, beta_y)

envelope_match_cache  = function_cache.Function_cache("envelope_match.cache")
def envelope_match(emit,current,g):
    if envelope_match_cache.in_cache(emit,current,g.get_mad_file(),
                                     g.get_line_name()):
        return envelope_match_cache.get(emit,current,g.get_mad_file(),
                                        g.get_line_name())
    (alpha_x, alpha_y, beta_x, beta_y) = get_alpha_beta(g)
    (s,kx,ky) = g.get_strengths()
    o = octapy.Octave()
    o.execute('LOADPATH="%s:";' % './envelope')
    o.set_value("alphax",alpha_x)
    o.set_value("alphay",alpha_y)
    o.set_value("betax",beta_x)
    o.set_value("betay",beta_y)
    o.set_value("s_array",s)
    o.set_value("Kx_array",kx)
    o.set_value("Ky_array",ky)
    o.set_value("emit",emit)
    o.set_value("current",current)
#    o.execute('printf("loadpath = %s\\n",LOADPATH)')
    o.execute('[sigma_x,sigma_xprime,r_x,sigma_y,sigma_yprime,r_y] = envelope_match(alphax,alphay,betax,betay,s_array,Kx_array,Ky_array, 0.4, current, emit, emit, 4, 1.0e-13)' )
    sigma_x = o.get_value("sigma_x")
    sigma_xprime = o.get_value("sigma_xprime")
    r_x = o.get_value("r_x")
    sigma_y = o.get_value("sigma_y")
    sigma_yprime = o.get_value("sigma_yprime")
    r_y = o.get_value("r_y")
    retval = [sigma_x,sigma_xprime,r_x,sigma_y,sigma_yprime,r_y]
    print "retval = ",retval
    if retval.count(None) == 0:
        envelope_match_cache.add(retval,emit,current,g.get_mad_file(),
                                 g.get_line_name())
    return retval

if ( __name__ == '__main__'):
    if len(sys.argv) > 1:
        map_order = int(sys.argv[1])
    else:
        map_order = 1
    if len(sys.argv) > 2:
        show_plot = int(sys.argv[2])
    else:
        show_plot = 1
    t0 = time.time()
    current = 0.5
    kinetic_energy = 0.4
    mass = physics_constants.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 200.0e6
    width_x = 0.004
    emittance = 2.77950448954e-06
    pipe_radius = 0.04
    num_particles = 100
    
    dpop = 3.0e-4
    z_frac = 0.1

    ee = error_eater.Error_eater()
    ee.start()
    g = gourmet.Gourmet("booster_classic.lat","bcel02",kinetic_energy,
                        scaling_frequency,order=map_order)
    g.insert_space_charge_markers(1)
    g.generate_linear_maps()
    t0 = time.time()

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=0)

    (alpha_x, alpha_y, beta_x, beta_y) = get_alpha_beta(g)

    envelope_match(1.0e-7,1.5,g)
    sys.exit(1)
    pz = bp.get_gamma() * bp.get_beta() * bp.mass_GeV

    (width_x,width_xprime,r_x) = matching.match_twiss_emittance(emittance,
                                                                alpha_x,
                                                                beta_x)
    bp.x_params(sigma = width_x, lam = width_xprime * pz)
    
    (width_y,width_yprime,r_y) = matching.match_twiss_emittance(emittance,
                                                                alpha_y,
                                                                beta_y)

    bp.y_params(sigma = width_y, lam = width_yprime * pz)

    sigma_z_meters = z_frac * bp.get_beta()*physics_constants.PH_MKS_c/scaling_frequency/math.pi
    bp.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    bp.correlation_coeffs(xpx = r_x, ypy = r_y)

    pgrid = processor_grid.Processor_grid(1)
    piperad = 0.04
    b = bunch.Bunch(current, bp, num_particles, pgrid)
    t0b = time.time()
    b.generate_particles()
    print "generate time =", time.time() - t0b
    

    the_map = g.get_linear_map(0)
    print Numeric.array2string(the_map,
                               precision=3,suppress_small=1)
    color = ['g','r','c','y','m','k']

    print g.get_u()
    mappers.crap(b.particles(),b.num_particles_local(), the_map)
    xsum = Numeric.zeros(num_particles,'d')
    xpsum = Numeric.zeros(num_particles,'d')
    ysum = Numeric.zeros(num_particles,'d')
    ypsum = Numeric.zeros(num_particles,'d')
    xnum = 0
    for count in range(0,100):
        g.get_fast_mapping(0).apply(b.particles(), b.num_particles_local())
#         mappers.apply_linear_map(b.particles(),b.num_particles_local(),
#                                  the_map)
        b.write_fort(count)
        xsum += b.particles()[0,:]
        xpsum += b.particles()[1,:]
        ysum += b.particles()[2,:]
        ypsum += b.particles()[3,:]
        xnum += 1
#         plot_bunch.plot_x(b,'b')
#         plot_bunch.plot_x(b,'g',1)
#         xsum0 += b.particles()[0,0]
#         xsum1 += b.particles()[0,1]
#         xnum += 1
            
    print "elapsed time =",time.time() - t0
    xav = xsum/xnum
    xpav = xpsum/xnum
    yav = ysum/xnum
    ypav = ysum/xnum
    xcalc = b.particles()[5,:]*the_map[0,5]
    xpcalc = b.particles()[5,:]*the_map[1,5]
    if show_plot:
        d = diagnostics.Diagnostics()
        pylab.plot(d.s,d.std[:,diagnostics.x],'r-')
        pylab.plot(d.s,d.std[:,diagnostics.y],'b-')
#        print "xav =",xsum/xnum, "xcalc =",b.particles()[5,:]*the_map[0,5]
        print the_map[0,5]
        print the_map[1,5]
        for i in range(0,10):
            print "%10.4g" % xav[i],
            print "%10.4g" % xpav[i],
            print "%10.4g" % yav[i],
            print "%10.4g" % ypav[i],
            print "%10.4g" % (xav[i]/xcalc[i]),
            print "%10.4g" % (xpav[i]/xpcalc[i])
#         print "xav1 =",xsum1/xnum, "xcalc1 =",b.particles()[5,1]*the_map[0,5]
        pylab.show()

