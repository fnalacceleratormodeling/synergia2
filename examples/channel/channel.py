#!/usr/bin/env bwpython

import Numeric
import time
import math
import os
import sys

import synergia
import s2_fish

from mpi4py import MPI

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
    solver = sys.argv[2]
    
    print "num_particles =",num_particles
    print "We will use a", solver, "solver"
    xwidth=0.0012026
    xpwidth=0.0049608
    rx=0.85440
    dpop = 1.0e-20

    ee = synergia.Error_eater()
    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.mad"),"channel",kinetic_energy,
                        scaling_frequency)
    gourmet.insert_space_charge_markers(kicks_per_line)

    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)
    print "Beam Beta", beam_parameters.get_beta()
    print "Beam Gamma", beam_parameters.get_gamma()
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 
    print "Betagamma and inverse betagamma",betagamma,1./betagamma
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = -rx)
    beam_parameters.y_params(sigma = xwidth, lam = xpwidth * pz,r = rx)
    sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
    beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)

    sys.stdout.flush()
    
    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)
    
    line_length = gourmet.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    diag = synergia.Diagnostics(gourmet.get_initial_u())
    kick_time = 0.0
    
    if solver == "3D" or solver == "3d":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_s2_fish=True,periodic=True)
    elif solver =="2D" or solver == "2d":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_gauss=True)
    print "elapsed time =",time.time() - t0

    diag.write_hdf5("channel")
    import pylab

    dimpact = synergia.Diagnostics_impact_orig("channel_impact_open")
    d0 = synergia.Diagnostics_impact_orig("channel0current")

    pylab.plot(d0.s,d0.std[:,synergia.x],'gx',label='no space charge')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    pylab.plot(dimpact.s,dimpact.std[:,synergia.x],'o',label='impact')
    pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],'r+',label='fish')
    pylab.legend(loc=0)
    pylab.show()
