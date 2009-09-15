#!/usr/bin/env bwpython

import numpy
import time
import math
import os
import sys

import synergia
import s2_fish
import dgourmet

from mpi4py import MPI
from pardebug import pardebug

if ( __name__ == '__main__'):
    t0 = time.time()

    myopts = synergia.Options("circular")
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",1,"map order",int)
    myopts.add("emittance95_over_pi",20,"95\% emittance in units of pi * mrad*mm",float)
    # longitudinal beta is 143.6
    myopts.add("dpop",3.482e-4,"(delta p)/p RMS width",float)
    #~ myopts.add("bunchlen", 0.05, "RMS bunchs length (z width) [m]", float)
    myopts.add("bunchlen", 40.0, "RMS bunch length (z width) [nanoseconds]", float)
    myopts.add("dpopoffset", 0.0, "offset in dpop", float)
    myopts.add("kicks",240,"kicks per line",int)
    myopts.add("turns",10,"number of turns",int)
    #~ myopts.add("latticefile","Debunch_modified.lat","",str)
    myopts.add("tgridnum",16,"transverse grid cells",int)
    myopts.add("lgridnum",16,"longitudinal grid cells",int)
    myopts.add("space_charge",1,"include space charge",int)
    myopts.add("impedance",0,"include impedance",int)
    myopts.add("energy",8.87710994,"total energy, default taken from Debunch_modified.lat",float)
    myopts.add("partpercell",1,"particles per cell",float)
    #~ myopts.add("current",0.5,"current",float)
    myopts.add("realnum",1.2e13,"number of real particles per bunch",float)
    myopts.add("solver","3d","solver",str)
    myopts.add("aperture",0.05,"aperture radius in m",float)
    myopts.add("numtrack",0,"number of particles to track",int)
    myopts.add("rampturns",30,"sextupole ramping turns",int)
    myopts.add("periodic",0,"longitudinal periodic boundary conditions",int)
    myopts.add("transverse",0,"use longitudinally uniform beam",int)
    myopts.add("tuneh",9.745,"horizontal fractional tune",float)
    myopts.add("tunev",9.93,"vertical fractional tune",float)
    
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      ["Debunch_modified.lat","dgourmet.py","debuncher.so"])

    t0 = time.time()
    if MPI.COMM_WORLD.Get_rank() ==0:
        log = open("log","w")
        output = "start"
        print output
        log.write("%s\n" % output)
        log.flush()
    energy = myopts.get("energy")
    mass = synergia.PH_NORM_mp
    kinetic_energy = energy-mass
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 2.35e6
    part_per_cell = myopts.get("partpercell")
    kicks_per_line = myopts.get("kicks")
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell)
    
 #   pipe_conduct= 1.4e6 # [ohm^-1 m^-1] (stainless steel)
    
    impedance=myopts.get("impedance")
    space_charge=myopts.get("space_charge")

    if MPI.COMM_WORLD.Get_rank() ==0:
        print "space_charge =",space_charge
        print "impedance =", impedance
        print "num_particles =",num_particles

    ee = synergia.Error_eater()
    ee.start()
    if MPI.COMM_WORLD.Get_rank() ==0:
        write_output = True
    else:
        write_output = False
            
    gourmet = dgourmet.Dgourmet(kinetic_energy,
                        scaling_frequency, myopts.get("maporder"),
                        myopts.get("tuneh"),myopts.get("tunev"),
                        write_output=write_output)
    
    #~ print numpy.array2string(gourmet.get_single_linear_map()[0:6,0:6],precision=2)
###    gourmet.print_elements(open("elements-orig.txt","w"))
    gourmet.insert_space_charge_markers(kicks_per_line)
    gourmet.complete_setup()
    orig_sextupoles = gourmet.get_sextupoles()
    if write_output:
        print "orig_sextupoles =",orig_sextupoles
    new_sextupoles = orig_sextupoles*0.0
    gourmet.set_sextupoles(new_sextupoles)
    if write_output:
        print "resulting sextupoles =",gourmet.get_sextupoles()
###    gourmet.print_elements(open("elements-withsc.txt","w"))
    #~ print numpy.array2string(gourmet.get_single_linear_map()[0:6,0:6],precision=2)
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    #~ print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
    
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=myopts.get("transverse"), adjust_zlength_to_freq=1)
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 

    emittance = myopts.get("emittance95_over_pi")*math.pi*1.0e-6/(6.0*math.pi*betagamma)
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV

    (xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emittance,alpha_x,beta_x)
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,
                             r = rx)
    if write_output:
        print "xwidth =",xwidth
    epsstar = emittance*4
    if write_output:
        print "eps* =",epsstar
    beta = beam_parameters.get_beta()
    if write_output:
        print "beta =",beta
    gamma = beam_parameters.get_gamma()
    if write_output:
        print "gamma =",gamma
    rp = 1.534698e-18
    if write_output:
        print "tune shift =",myopts.get("realnum")*rp/(math.pi * epsstar * beta *gamma *gamma)
    
    (ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emittance,alpha_y,beta_y)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,
                             r = ry)
    
    bunchlen_sec = myopts.get("bunchlen")*1e-9
    bunchlen_m = bunchlen_sec * beam_parameters.get_beta() * synergia.PH_MKS_c
    if write_output:
        print "bunchlen_m =", bunchlen_m
    beam_parameters.z_params(sigma = bunchlen_m,
                             lam = myopts.get("dpop")* pz)
    
    sys.stdout.flush()
    
    s = 0.0
    line_length = gourmet.orbit_length()
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "line_length =",line_length
    if kicks_per_line > 0:
        tau = 0.5*line_length/kicks_per_line
    else:
        tau = 1.0 # not used!!!
    kick_time = 0.0
    beta = beam_parameters.get_beta()
    
    bunchnp = 1.0e9;
    bunch = s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,periodic=True)
    
    #current = myopts.get("realnum")*synergia.PH_MKS_e*scaling_frequency
    #if write_output:
    #    print "current =",current
    
    #bunch = s2_fish.Macro_bunch(mass,1)
    #bunch.init_gaussian(num_particles,current,beam_parameters)
    bunch.write_particles("begin")

    if myopts.get("numtrack") > 0:
        tracker = synergia.Tracker("/tmp",(myopts.get("numtrack"),num_particles))
        tracker.add(bunch,0.0)
    else:
        tracker = None
    track_period_steps=80

    diag = synergia.Diagnostics(gourmet.get_initial_u(),short=True)

    if MPI.COMM_WORLD.Get_rank() ==0:
        output = "setup time = %g" % (time.time() - t0)
        print output
        log.write("%s\n" % output)
        log.flush()
    t5total = 0
    
    if myopts.get("solver") == "3d" or myopts.get("solver") == "3D":
        space_charge_solver ="s2_fish_3d"
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "using 3d solver"
    elif myopts.get("solver") == "2d" or myopts.get("solver") == "2D":
        space_charge_solver = "s2_fish_transverse"
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "using 2d transverse solver"
    else:
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "unknown solver",myopts.get("solver")
        sys.exit(1)
    
    space_charge = s2_fish.SpaceCharge(space_charge_solver, grid=griddim, periodic=myopts.get("periodic"), transverse=True)    
        
    for turn in range(1,myopts.get("turns")+1):
        t1 = time.time()
        if turn<= myopts.get("rampturns"):
            new_sextupoles = orig_sextupoles*turn/(1.0*myopts.get("rampturns"))
            gourmet.set_sextupoles(new_sextupoles)
            gourmet.generate_actions()

        s = synergia.propagate(s,gourmet,bunch,diag, 
                                impedance=impedance, space_charge=space_charge,aperture=myopts.get("aperture"), 
                                tracker=tracker,track_period_steps=track_period_steps)
        ### keep beam from moving
#         local_sump = numpy.zeros([6],'d')
#         global_sump = numpy.zeros([6],'d')
#         for i in range(0,6):
#             local_sump[i] = numpy.sum(bunch.get_local_particles()[i,:])
#         MPI.COMM_WORLD.Allreduce([local_sump,MPI.DOUBLE],
#                                  [global_sump,MPI.DOUBLE],op=MPI.SUM)
#         newmeansp = global_sump/bunch.total_num

#         for i in range(0,6):
#             bunch.get_local_particles()[i,:] += -newmeansp[i]
#         for i in range(0,6):
#             local_sump[i] = numpy.sum(bunch.get_local_particles()[i,:])
#         MPI.COMM_WORLD.Allreduce([local_sump,MPI.DOUBLE],
#                                  [global_sump,MPI.DOUBLE],op=MPI.SUM)
#         finalmeansp = global_sump/bunch.total_num
#         pardebug("new: %s\nfinal:%s\n" % (newmeansp,finalmeansp))
        ### end keep beam from moving
        t2 = time.time()
        #~ mbs = bunch.get_store()
        #~ s2_fish.constraints.apply_circular_aperture(mbs,myopts.get("aperture"))
        #~ # egad
        #~ bunch.local_num = mbs.local_num
        #~ bunch.total_num = mbs.total_num
        t3 = time.time()
        bunch.write_particles("turn_%03d.h5" % turn)
        t4 = time.time()
        if MPI.COMM_WORLD.Get_rank() ==0:
            output = "turn %d time = %g, %g %g %g %g" % (turn, t4 -t1, t2 - t1, t3-t2, t4-t3, t5total)
            # times: 0id, 1allbutlog 2propagate 3constraint 4write 5log
            print output
            log.write("%s\n" % output)
            log.flush()
        t5total = time.time() - t4 
    if MPI.COMM_WORLD.Get_rank() == 0:
        diag.write_hdf5("diagnostics")
    bunch.write_particles("end")
    if tracker:
        tracker.close()
        tracker.show_statistics()    
    if MPI.COMM_WORLD.Get_rank() ==0:
        log.close()
        print "elapsed time =",time.time() - t0
 
 
