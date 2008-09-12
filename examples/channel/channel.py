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
    if len(sys.argv) <=2:
        print "usage: channel.py gridnum solver [xoffset] [impedance] [pipe_radius] [pipe_conduct] [space_charge]"
        sys.exit(1)

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
    xoffset = 0.0
    if len(sys.argv)>3:
        xoffset = float(sys.argv[3])
    
    impedance = False
    if len(sys.argv)>4:
        if (sys.argv[4] == "True" or sys.argv[4] == "true"):
            impedance = True
        elif (sys.argv[4] == "False" or sys.argv[4] == "false"):
            impedance = False
        else:
            raise RuntimeError,\
                "impedance (third argument) must be true or false"
    
    pipe_radius = 0.01
    if len(sys.argv)>5:
        pipe_radius = float(sys.argv[5])
        
    pipe_conduct= 1.4e6 # [/s] (stainless steel)
    if len(sys.argv)>6:
        pipe_conduct= float(sys.argv[6])
    
    space_charge = True
    if len(sys.argv)>7:
        if (sys.argv[7] == "True" or sys.argv[7] == "true"):
            space_charge = True
        elif (sys.argv[7] == "False" or sys.argv[7] == "false"):
            space_charge = False
        else:
            raise RuntimeError,\
                "space_charge (seventh argument) must be true or false"
    
    print "space_charge =",space_charge
    
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
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = -rx,offset=xoffset)
    beam_parameters.y_params(sigma = xwidth, lam = xpwidth * pz,r = rx)
    sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
    beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)

    sys.stdout.flush()
    
    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)
    bunch.write_particles("begin")
    line_length = gourmet.orbit_length()
    tau = 0.5*line_length/kicks_per_line
    s = 0.0
    first_action = 0
    diag = synergia.Diagnostics(gourmet.get_initial_u())
    kick_time = 0.0
    
    if solver == "3D" or solver == "3d":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_s2_fish=True,periodic=True,
            impedance=impedance,space_charge=space_charge,
            pipe_radiusx=pipe_radius,pipe_radiusy=pipe_radius, pipe_conduct=pipe_conduct)
    if solver == "3DC" or solver == "3dc":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,
            use_s2_fish_cylindrical=True,radius=0.01,
            impedance=impedance,space_charge=space_charge,
            pipe_radius=pipe_radius,pipe_conduct=pipe_conduct)
    elif solver =="2D" or solver == "2d":
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_gauss=True,
            impedance=impedance,pipe_radius=pipe_radius,pipe_conduct=pipe_conduct)
    print "elapsed time =",time.time() - t0
    bunch.write_particles("end")
    diag.write_hdf5("channel")
    import pylab

    dimpact = synergia.Diagnostics_impact_orig("channel_impact_open")
    d0 = synergia.Diagnostics_impact_orig("channel0current")

    pylab.plot(d0.s,d0.std[:,synergia.x],'gx',label='no space charge')
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')

    pylab.plot(dimpact.s,dimpact.std[:,synergia.x],'o',label='impact')
    pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],'r+',markersize=15.0,label='fish')
    pylab.legend(loc=0)
    pylab.show()
