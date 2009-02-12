#!/usr/bin/env bwpython

import local_paths
import gourmet
import numpy
import numpy
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

import pylab

#import do_compare_channel
import memory

import UberPkgpy #SpaceChargePkgpy

import mappers
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

def plot_bunch(b,color):
    beg1 = 0
    end1 = 1000
    pylab.subplot(2,2,1)
    pylab.title('x')
    x1 = b.particles()[0,beg1:end1]
    xprime1 = b.particles()[1,beg1:end1]

    end2 = len(b.particles()[0,:])
    beg2 = end2 - (end1 - beg1)
    x2 = b.particles()[0,beg2:end2]
    xprime2 = b.particles()[1,beg2:end2]
    
    pylab.plot(x1,xprime1,'%s.' % color)
    pylab.plot(x2,xprime2,'%s.' % color)
    
    pylab.subplot(2,2,2)
    pylab.title('y')
    y1 = b.particles()[2,beg1:end1]
    yprime1 = b.particles()[3,beg1:end1]

    y2 = b.particles()[2,beg2:end2]
    yprime2 = b.particles()[3,beg2:end2]

    pylab.plot(y1,yprime1,'%s.' % color)
    pylab.plot(y2,yprime2,'%s.' % color)

    pylab.subplot(2,2,3)
    pylab.title('z')
    z1 = b.particles()[4,beg1:end1]
    zprime1 = b.particles()[5,beg1:end1]

    z2 = b.particles()[4,beg2:end2]
    zprime2 = b.particles()[5,beg2:end2]
    

    pylab.plot(z1,zprime1,'%s.' % color)
    pylab.plot(z2,zprime2,'%s.' % color)

mt = 0
def wrapped_apply_map(map, the_bunch):
    global mt
    t0 = time.time()
    apply_map.apply_map1(the_bunch.particles(),the_bunch.num_particles_local(),
                         map)
    mt += time.time() - t0

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
    num_particles = 1000000
    
    dpop = 3.0e-4
    z_frac = 0.1

    ee = error_eater.Error_eater()
    ee.start()
    g = gourmet.Gourmet("simplebooster.mad","fullb",kinetic_energy,
                        scaling_frequency,order=map_order)
    g.insert_space_charge_markers(1)
    g.generate_linear_maps()
    t0 = time.time()
    print "chef map generation took", time.time()-t0

    bp = beam_parameters.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=0)

    (beta_x, beta_y, alpha_x, alpha_y) =  (5.756421714817864,
                                           20.29327699231588,
                                           5.209874098706389e-15,
                                           5.796656894202063e-15)
    
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
    bp.correlation_coeffs(xpx = r_x, ypy = r_x)

    pgrid = processor_grid.Processor_grid(1)
    piperad = 0.04
    b = bunch.Bunch(current, bp, num_particles, pgrid)
    t0b = time.time()
    b.generate_particles()
    print "generate time =", time.time() - t0b
    
    plot_bunch(b,'b')

    color = ['g','r','c','y','m','k']
    map_times = []
    convert_times = []
    print g.get_linear_map(0)
    for count in range(0,2):
###        g.get_fast_mapping(0).apply(b.particles(), b.num_particles_local())
        mappers.apply_linear_map(b.particles(),b.num_particles_local(),
                                 g.get_linear_map(0))
        plot_bunch(b,color[count])
            
    print "elapsed time =",time.time() - t0
    if show_plot:
        pylab.show()

