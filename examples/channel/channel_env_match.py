#!/usr/bin/env bwpython

import numpy
import time
import math
import os
import sys

import synergia
import s2_fish

from mpi4py import MPI

if ( __name__ == '__main__'):
    t0 = time.time()
    myopts = synergia.Options("channel")
    myopts.add("current",0.5,"current",float)
    myopts.add("kinetic_energy",0.0067,"kinetic energy in GeV",float)
    myopts.add("Solver","2D","Solver 3D or 2D",str)
    myopts.add("griddim",[16,16,16],"Space charge grid for 3D",int)
    myopts.add("part_per_cell",1,"particlels per cell",int)
    
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)

    mass = synergia.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6
    pipe_radius = 0.04
    kicks_per_line = 10

    griddim=myopts.get("griddim")
    part_per_cell=myopts.get("part_per_cell")
    num_particles = int(griddim[0]*griddim[1]*griddim[2]) * part_per_cell
    
    print "num_particles =",num_particles
    Solver=myopts.get("Solver")
    print "We will use a ", Solver, " solver"

    # these are the values for a matched 0.0067, 0.5 channel
    xwidth=0.0012026
    xpwidth=0.0049608
    rx=0.85440
    # end of block

    kinetic_energy = myopts.get("kinetic_energy")
    current = myopts.get("current")
    dpop = 1.0e-20

    ee = synergia.Error_eater()
    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.mad"),
                               "channel",kinetic_energy,scaling_frequency)
    gourmet.insert_space_charge_markers(kicks_per_line)

    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)
    print "Beam Beta", beam_parameters.get_beta()
    print "Beam Gamma", beam_parameters.get_gamma()
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 
    print "Betagamma and inverse betagamma",betagamma,1./betagamma

    #matching attempt
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.get_alpha_beta(gourmet)
    #Now use the x_width we had for the low energy case, to get emittance 
    (xpwidth,rx,emittancex) = synergia.match_twiss_width(xwidth,alpha_x,beta_x)
    (ypwidth,ry,emittancey) = synergia.match_twiss_width(xwidth,alpha_y,beta_y)
    (sigma_x,sigma_xprime,rx,sigma_y,sigma_yprime,ry) = \
                                                      synergia.envelope_match(emittancex,emittancey,current,gourmet)
    print "sigma_x,sigma_xprime,rx,sigma_y,sigma_yprime,ry ", sigma_x,sigma_xprime,rx,sigma_y,sigma_yprime,ry
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = sigma_x, lam = sigma_xprime * pz,r = -rx)
    beam_parameters.y_params(sigma = sigma_y, lam = sigma_yprime * pz,r = -ry)
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
    
    if Solver == "3D":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_s2_fish=True)
    elif Solver =="2D":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_gauss=True)
    print "elapsed time =",time.time() - t0

    diag.write_hdf5("channel")
    import pylab

    dimpact = synergia.Diagnostics_impact_orig("channel_impact_open")
    d0 = synergia.Diagnostics_impact_orig("channel0current")

    pylab.plot(d0.s,d0.std[:,synergia.x],'gx',label='no SC, 0.0067 GeV')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    pylab.plot(dimpact.s,dimpact.std[:,synergia.x],'o',label='impact, 0.0067 GeV')
    pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],'r+',label='fish')
    pylab.legend(loc=0)
    envelope = synergia.loadfile("envelope_match.dat")
    pylab.plot(envelope[:,0],envelope[:,1])
    pylab.show()
