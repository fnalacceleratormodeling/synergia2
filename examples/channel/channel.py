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
    myopts.add("gridnum",16,"number of grid points to be used for all directions",int)
    myopts.add("solver","3d","solver",str)
    myopts.add("xoffset",0.0,"x offset",float)
    myopts.add("impedance",0,"whether to use resistive wall kicks",int)
    myopts.add("piperadius",0.01,"pipe radius for impedance",float)
    myopts.add("pipeconduct",1.4e6,
        "conductivity for pipe [/s], default is for stainless steel",float)
    myopts.add("spacecharge",1,"whether to use space charge kicks",int)        
    myopts.add("doplot",1,"show plot",int)
    myopts.add("current",0.5,"current -- modify at your own risk",float)

    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      ["channel.mad"])    
    
    current = myopts.get("current")
    kinetic_energy = 0.0067
    mass = synergia.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6
    part_per_cell = 1
    width_x = 0.004
    kicks_per_line = 30
    gridnum = myopts.get("gridnum")
    griddim = (gridnum,gridnum,gridnum)
    num_particles = griddim[0]*griddim[1]*griddim[2] * part_per_cell

    xoffset = myopts.get("xoffset")  
    pipe_radius = myopts.get("piperadius")
    pipe_conduct= myopts.get("pipeconduct")
    space_charge = myopts.get("spacecharge")
    solver = myopts.get("solver")
    impedance = myopts.get("impedance")
    
    xwidth_initial=0.0012026
    ywidth_initial=0.0012026
    #xwidth=0.0012026
    #xpwidth=0.0049608
    #rx=0.85440
    dpop = 1.0e-20

    ee = synergia.Error_eater()
    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.xsif"),"channel",kinetic_energy,
                        scaling_frequency)
    gourmet.insert_space_charge_markers(kicks_per_line) 
    
  
    
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 

    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    #~print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
    emitx=(xwidth_initial)**2/beta_x
    emity=(ywidth_initial)**2/beta_y
    print "emitx= ", emitx, "emity = ", emity
   
    
  
    xwidth=xwidth_initial
    xpwidth=0.0049608
    rx=-0.85440
   
   
    ywidth=ywidth_initial
    ypwidth=0.0049608
    ry=0.85440
    
##~~   input twiss matched beam parameters, require only the emitance   
    #(xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emitx,alpha_x,beta_x)
    #(ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emity,alpha_y,beta_y)
     #~retval=synergia.matching.envelope_match(emitx,emity,current,gourmet,do_plot=1) 
    
    
    widths=[xwidth,xpwidth,rx,ywidth,ypwidth,ry]
    retval=synergia.matching.envelope_motion(widths,current,gourmet,do_plot=1,do_match=0)
   
    
    #betax_chef=gourmet.get_lattice_functions().beta_x 
    #s_chef =gourmet.get_lattice_functions().s
    
    #sigmax_chef=numpy.zeros(len(betax_chef),'d')
    #for i in range(len(betax_chef)):
	    #sigmax_chef[i]=math.sqrt(betax_chef[i]*emitx)
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "num_particles =",num_particles
        print "We will use a", solver, "solver"
        print "Beam Beta", beam_parameters.get_beta()
        print "Beam Gamma", beam_parameters.get_gamma()
        print "Betagamma and inverse betagamma",betagamma,1./betagamma
    
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry)
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
    
    use_s2_fish = False
    use_s2_fish_cylindrical = False
    use_gauss = False
    radius=None
    if solver == "3D" or solver == "3d":
        use_s2_fish = True
    if solver == "3DC" or solver == "3dc":
        use_s2_fish_cylindrical = True
        radius = 0.01
    elif solver =="2D" or solver == "2d":
        use_gauss = True
    s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,
                           use_s2_fish=use_s2_fish,
                           use_s2_fish_cylindrical=use_s2_fish_cylindrical,
                           use_gauss=use_gauss,
                           radius=radius,
                           periodic=True,
                           impedance=impedance,space_charge=space_charge,
                           pipe_radiusx=pipe_radius,pipe_radiusy=pipe_radius,
                           pipe_conduct=pipe_conduct)
    print "elapsed time =",time.time() - t0,"on rank", MPI.COMM_WORLD.Get_rank()
    bunch.write_particles("end")
    diag.write_hdf5("channel")
    if myopts.get("doplot") and MPI.COMM_WORLD.Get_rank() == 0:
        import pylab
	#pylab.plot(s_chef, sigmax_chef,'yp',label='chef no sp ch')   
	
	
	t0=time.time()
	solver="2d"
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_gauss=True,
            impedance=impedance,pipe_radiusx=pipe_radius,pipe_radiusy=pipe_radius,pipe_conduct=pipe_conduct)
    
        print "elapsed time 2d =",time.time() - t0,"on rank", MPI.COMM_WORLD.Get_rank()
              
	
	dimpact = synergia.Diagnostics_impact_orig("channel_impact_open")
        d0 = synergia.Diagnostics_impact_orig("channel0current")

	pylab.figure(1)
        pylab.plot(d0.s,d0.std[:,synergia.x],'gx',label='no space charge')
        pylab.xlabel('s (m)')
        pylab.ylabel('std<x> (m)')
        pylab.plot(dimpact.s,dimpact.std[:,synergia.x],'o',label='impact')
        pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],'r+',markersize=15.0,label='fish 2d')
        
	
	pylab.figure(2)
	pylab.plot(d0.s,d0.std[:,synergia.xprime],'gx',label='no space charge')
	pylab.plot(dimpact.s,dimpact.std[:,synergia.xprime],'o',label='impact')
	pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.xprime],'r+',markersize=15.0,label='fish 2d')	
        pylab.xlabel('s (m)')
        pylab.ylabel('std<xprime> ')
	
	
	solver="3d"
	
	bunch = s2_fish.Macro_bunch(mass,1)
        bunch.init_gaussian(num_particles,current,beam_parameters)
        bunch.write_particles("begin")
        s = 0.0
        first_action = 0
        diag = synergia.Diagnostics(gourmet.get_initial_u())
        kick_time = 0.0
    
	
	pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
        beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
        beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry)
        sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
        beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)
	
	t0=time.time()
        s = synergia.propagate(0.0,gourmet,bunch,diag,griddim,use_s2_fish=True,periodic=True,
            impedance=impedance,space_charge=space_charge,
            pipe_radiusx=pipe_radius,pipe_radiusy=pipe_radius, pipe_conduct=pipe_conduct)
        print "elapsed time 3d =",time.time() - t0,"on rank", MPI.COMM_WORLD.Get_rank()
	bunch.write_particles("end")
        diag.write_hdf5("channel")
	
	pylab.figure(1)
	pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],'b+',label='fish 3d')
	pylab.legend(loc=0)
	
	pylab.figure(2)
	pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.xprime],'b+',label='fish 3d')
	pylab.legend(loc=0)
	
	pylab.show()
