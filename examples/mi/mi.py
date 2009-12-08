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

    myopts = synergia.Options("mi")
    #~ myopts.add("current",0.5,"current",float)
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("emittance",3.0e-6,"emittance",float)
    myopts.add("dpop",1.0e-20,"(delta p)/p",float)
    myopts.add("dpopoffset", 0.0, "offset in dpop (am: ! -duop)", float)
    myopts.add("kicks",10,"kicks per line",int)
    myopts.add("turns",4,"number of turns",int)
    myopts.add("latticefile","mi20-egs.lat","",str)
    myopts.add("tgridnum",16,"transverse grid cells",int)
    myopts.add("lgridnum",64,"",int)
    myopts.add("xoffset",0,"",float)
    myopts.add("yoffset",0,"",float)
    myopts.add("zoffset",0.1,"offset in z", float)
    myopts.add("xpoffset", 0, "offset in x-prime", float)
    myopts.add("ypoffset", 0, "offset in y-prime", float)
    myopts.add("emittance",0.26e-6,"",float)
    myopts.add("space_charge",1,"",int)
    myopts.add("impedance",0,"",int)
    myopts.add("energy",8.9,"",float)
    myopts.add("partpercell",1,"",float)
    myopts.add("bunches",1,"",int)
    myopts.add("bunchnp",1.0e11,"number of particles per bunch",float)
    myopts.add("pipe_radius", 0.025, "pipe radius", float)
    
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      [myopts.get("latticefile")])

    t0 = time.time()
    energy = myopts.get("energy")
    mass = synergia.PH_NORM_mp
    kinetic_energy = energy-mass
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 53e6
    pipexradius = 0.123
    pipeyradius = 0.0508
    part_per_cell = myopts.get("partpercell")
    kicks_per_line = myopts.get("kicks")
    xoffset = myopts.get("xoffset")
    yoffset = myopts.get("yoffset")

    pipe_conduct= 1.4e6 # [/s] (stainless steel)
    
    impedance=myopts.get("impedance")
    space_charge=myopts.get("space_charge")


    ee = synergia.Error_eater()
    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile"))
        ,"ring_p_q605",kinetic_energy,
                        scaling_frequency)
    gourmet.insert_space_charge_markers(kicks_per_line)
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    #~ print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
    
    
    # defining beam_parameters
    charge = 1.0
    initial_phase = 0.0

    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency_Hz=scaling_frequency,
                                         transverse=0, adjust_zlength_to_freq=0)

    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 

    emittance = myopts.get("emittance")
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV

    (xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emittance,alpha_x,beta_x)
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
    
    (ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emittance,alpha_y,beta_y)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry,offset=yoffset)
    
    sigma_z_meters = 0.75
    beam_parameters.z_params(sigma = sigma_z_meters, lam = myopts.get("dpop")* pz)

    sys.stdout.flush()
    
    s = 0.0
    line_length = gourmet.orbit_length()
    bunch_spacing = line_length/588.0

    if MPI.COMM_WORLD.Get_rank() ==0:
        print "space_charge =",space_charge
        print "impedance =",impedance
        
        print "line_length =",line_length
        print "bunch_spacing =",bunch_spacing


    tau = 0.5*line_length/kicks_per_line
    kick_time = 0.0
    beta = beam_parameters.get_beta()
    bunchnp = myopts.get("bunchnp")
    current = bunchnp * \
        synergia.physics_constants.PH_MKS_e/ \
        (bunch_spacing/(beta*synergia.physics_constants.PH_MKS_c))
    print "current =",current, " (but not used)"
    
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell)

    # creating the bunch
    numbunches = myopts.get("bunches") # I can only deal with one bunch now
    bunches = []
    for bunchnum in range(0,numbunches):
        diag=synergia.Diagnostics(gourmet.get_initial_u(),short=True)
        bunches.append(s2_fish.Macro_bunch.gaussian(bunchnp, num_particles,
                                                    beam_parameters,
                                                    bucket_num = 0,
                                                    diagnostics=diag,
                                                    periodic=False))
        bunches[bunchnum].write_particles("begin-%02d"%bunchnum)

    if space_charge:
        solver = "s2_fish_3d"
        sp_ch_obj =     beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency_Hz=scaling_frequency,
                                         transverse=0, adjust_zlength_to_freq=0)
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 

    emittance = myopts.get("emittance")
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV


    (xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emittance,alpha_x,beta_x)
    xoffset = myopts.get("xoffset")
    xpoffset = myopts.get("xpoffset")
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,
                             r = rx,offset=xoffset, offset_p = xpoffset * pz)
    
    
    yoffset = myopts.get("yoffset")
    ypoffset = myopts.get("ypoffset")
    (ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emittance,alpha_y,beta_y)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,
                             r = ry,offset=yoffset, offset_p = ypoffset * pz )
    
    
    #zoffset = myopts.get("zoffset")
    #zwidth=myopts.get("bunchlen")
    zpwidth=myopts.get("dpop")
    
    bunch_sp=2.0*math.pi*beta*synergia.physics_constants.PH_MKS_c/beam_parameters.get_omega()
    z_length=bunch_sp
    zwidth=bunch_sp/15.0    
    zoffset = 0.2
    beam_parameters.z_params(sigma = zwidth, 
                             lam = zpwidth* pz, z_length=z_length, offset=zoffset,
                             offset_p = myopts.get("dpopoffset")*pz)
#Note! the input term  dpop is in fact (-duop)                           
  

### creating the bunch
    bunchnp0=myopts.get("bunchnp") 
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    part_per_cell = myopts.get("partpercell")
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell)
    
    numbunches = myopts.get("bunches")
    bunches = []

   
    if MPI.COMM_WORLD.Get_rank() ==0:
        print " **********************************************************************" 
        print "beam information:"
        print
        print "transverse=",beam_parameters.transverse
        if beam_parameters.adjust_zlength_to_freq:
            print "adjusted length =",beam_parameters.z_length
        else:
            print "bunch length (not adjusted)= ", beam_parameters.z_length
           
        print "grid size: ", griddim
        print "number of macroparticles (num_particles): ", num_particles
        print "Charge of bunch: ", bunchnp0
        print "Number of bunches: ", numbunches
        
    sys.stdout.flush()

    for bunchnum in range(0,numbunches):
        diag=synergia.Diagnostics(gourmet.get_initial_u(),short=True)
        bunchnp=bunchnp0#*(bunchnum+1)*0.5 # bucket_num =2 in front of bucket_num =3
        bunches.append(s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,diagnostics=diag,bucket_num=2*bunchnum,periodic=True))
        bunches[bunchnum].write_particles("begin-%02d"%bunchnum)
        print " bunch(",bunchnum,") periodicity=",bunches[bunchnum].periodic
       # print "  initial means bunch(",bunchnum,")=",numpy.array(bunches[bunchnum].diagnostics.get_means())

    print " **********************************************************************"  
    mbunches=s2_fish.Multiple_bunches(bunches, bunch_sp)
        
    #diagnostic_units=gourmet.get_initial_u()
    #bunch = s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,gourmet.get_initial_u(),periodic=True)
    #bunch= s2_fish.Macro_bunch.from_bunch(bunch1)
    #bunch = s2_fish.Macro_bunch.test(int(part_per_cell))
   # bunch = s2_fish.Macro_bunch.test_am(bunchnp,part_per_cell,griddim,beam_parameters)
   # bunch = s2_fish.Macro_bunch.sphere(num_particles, 0.001)
   # bunch = s2_fish.Macro_bunch.cylinder(num_particles, 0.001,0.01)
   # covar=beam_parameters.get_covariances()
   # bunch = s2_fish.Macro_bunch.gaussian_covariance(bunchnp,num_particles,beam_parameters,covar,periodic=True)
   # print "bunch first long size=",bunch.get_longitudinal_period_size()
   
   
    #bunch.write_particles("begin")
    #diag = synergia.Diagnostics(gourmet.get_initial_u(),save_period=0)   

    log = open("log","w")
    if MPI.COMM_WORLD.Get_rank() ==0:
            output = "start propagation"
            print output
            log.write("%s\n" % output)
            log.flush()

 
    # solver and geometry

    pipe_radius = myopts.get("pipe_radius")
    space_charge=myopts.get("space_charge")
    if space_charge:
        solver="s2_fish_3d"
        sp_ch_obj = s2_fish.SpaceCharge(solver,griddim,radius_cylindrical=pipe_radius,periodic=True)	
        print " sp_ch grid=",sp_ch_obj.get_grid()
        print " sp_ch solver=",sp_ch_obj.get_solver()
        print " sp_ch pipe radius=",sp_ch_obj.get_radius_cylindrical()
    else:
	   sp_ch_obj=None    	

        
    log = open("log","w")
    if MPI.COMM_WORLD.Get_rank() ==0:
            output = "start propagation"
            print output
            log.write("%s\n" % output)
            log.flush()
    for turn in range(1,myopts.get("turns")+1):
        t1 = time.time()
        s = synergia.propagate(s,gourmet,
            bunches, space_charge=sp_ch_obj,
            impedance=None)
        if MPI.COMM_WORLD.Get_rank() ==0:
            output = "turn %d time = %g"%(turn,time.time() - t1)
            print output
            log.write("%s\n" % output)
            log.flush()
    for bunchnum in range(0,numbunches):
        if MPI.COMM_WORLD.Get_rank() == 0:
            bunches[bunchnum].diagnostics.write_hdf5("mi-%02d"%bunchnum)
    for bunchnum in range(0,numbunches):
        bunches[bunchnum].write_particles("end-%02d"%bunchnum)
    log.close()
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "elapsed time =",time.time() - t0
 
 
