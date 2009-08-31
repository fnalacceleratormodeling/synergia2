#!/usr/bin/env bwpython
import numpy
import time
import math
import os
import sys

import synergia
import s2_fish
import impact
import pylab

from mpi4py import MPI

SOLVER_ORBIT = "orbit"
SOLVER_GAUSS  = "gauss"
SOLVER_3D = "solver3d"
NO_SOLVER = "no_solver"

def runSolverTest(beam_parameters,mass,num_particles,current_in,gourmet, solver_type,color='k', periodic = True):
    current=current_in
    print "\n***************** testing ",solver_type
    print "num_particles =",num_particles
    print "griddim =",griddim
    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)
    #modify_test_particles(bunch)
    diag = synergia.Diagnostics(gourmet.get_initial_u())
    t0=time.time()
    
    if solver_type == NO_SOLVER:
        space_charge = None
        space_charge_solver = ""
    elif solver_type == SOLVER_3D:
        space_charge_solver="s2_fish_3d"   
    elif solver_type == SOLVER_ORBIT:
        space_charge_solver="s2_fish_transverse" 		
    elif solver_type == SOLVER_GAUSS:
        space_charge_solver="s2_fish_gauss"     
    else:
        print "ERROR: unknown solver!"
        return; 
    
    if space_charge_solver != "":
        space_charge = s2_fish.SpaceCharge(space_charge_solver, grid=griddim, periodic=periodic)

    s = synergia.propagate(0.0,gourmet, bunch,diag, space_charge=space_charge)
	  
    print "elapsed time =",time.time() - t0,"on rank", MPI.COMM_WORLD.Get_rank()
    solver_type_label=solver_type
    if current_in==0:
        solver_type_label="no spc"     
    pylab.figure(1)
    pylab.xlabel('s (m)')
    pylab.ylabel('std<x> (m)')	
    pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.x],color, label=solver_type_label)
    pylab.legend(loc=0)
    pylab.figure(2)
    pylab.xlabel('s (m)')
    pylab.ylabel('std<xprime> ')
    pylab.plot(diag.get_s(),diag.get_stds()[:,synergia.xprime],color, label=solver_type_label)
    pylab.legend(loc=0)
    return diag
  

# call this after init_gaussian
def modify_test_particles(bunch):
    num = 20
    max = 1.2
    means,std = synergia.get_spatial_means_stds(bunch)
    print "std ",std
    for i in range(0,num):
	bunch.get_local_particles()[0,i] = max*i/(1.0*num)
	bunch.get_local_particles()[1,i] = 0.0
	bunch.get_local_particles()[2,i] = max*0.5
	bunch.get_local_particles()[3,i] = 0.0

    
if ( __name__ == '__main__'):
    t0 = time.time()
    myopts = synergia.Options("channel")
    myopts.add("gridnum",32,"number of grid points to be used for all directions",int)
    myopts.add("solver","3d","solver",str)
    myopts.add("xoffset",0.0,"x offset",float)
    myopts.add("impedance",0,"whether to use resistive wall kicks",int)
    myopts.add("piperadius",100.,"pipe radius for impedance",float)
    myopts.add("pipeconduct",1.4e6,
        "conductivity for pipe [/s], default is for stainless steel",float)
    myopts.add("spacecharge",1,"whether to use space charge kicks",int)        
    myopts.add("kinetic",0.0067,"kinetic energy", float)
    myopts.add("current",0.1,"in current value",float)
    myopts.add("do_plot",1,"draw plot", int)
    myopts.add("x_initial",0.0012026,"initial width", float)
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,["channel.mad"])    
    
    current_in = myopts.get("current")
    do_plot = myopts.get("do_plot")
     
    print "curent=",current_in
    
    kinetic_energy = myopts.get("kinetic");
    print "kinetic_energy= ",kinetic_energy
    mass = synergia.PH_NORM_mp
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = 10221.05558e6
#    scaling_frequency = 47713451.5923694e3
    part_per_cell = 1
    kicks_per_line = 40
    gridnum = myopts.get("gridnum")

    xoffset = myopts.get("xoffset")  
    pipe_radius = myopts.get("piperadius")
    pipe_conduct= myopts.get("pipeconduct")
    space_charge = myopts.get("spacecharge")
    solver = myopts.get("solver")
    impedance = myopts.get("impedance")
    xwidth_initial = myopts.get("x_initial")
    #xwidth_initial=0.0012026
    ywidth_initial = xwidth_initial
    #ywidth_initial=0.0012026
    dpop = 1.0e-20

    ee = synergia.Error_eater()
    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),"channel.mad"),"channel",kinetic_energy,
                        scaling_frequency,delay_complete=True)
			
    gourmet.insert_space_charge_markers(kicks_per_line) 
    gourmet.complete_setup()
     
    print "line_length =", gourmet.orbit_length()
    
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=1)
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 
    beam_parameters.z_peaks = 1
    
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    emitx=(xwidth_initial)**2/beta_x
    emity=(ywidth_initial)**2/beta_y
    print "emitx= ", emitx, "emity = ", emity
  
    xwidth=xwidth_initial
    xpwidth=0.0049608
    rx=-0.85440
    
    ywidth=ywidth_initial
    ypwidth=0.0049608
    ry=0.85440
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "Beam Beta", beam_parameters.get_beta()
        print "Beam Gamma", beam_parameters.get_gamma()
        print "Betagamma and inverse betagamma",betagamma,1./betagamma
    
    radius=10.0*math.sqrt(xwidth*xwidth+ywidth*ywidth)
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry)
    sigma_z_meters = beam_parameters.get_beta()*synergia.PH_MKS_c/scaling_frequency/math.pi
    beam_parameters.z_params(sigma = sigma_z_meters, lam = dpop* pz)
    sys.stdout.flush()
  
    widths=[xwidth,xpwidth,rx,ywidth,ypwidth,ry]
    current=current_in
    synergia.matching.envelope_motion(widths,current,gourmet,do_plot=1,do_match=0)
   
    # Orbit solver    
    griddim = (gridnum,gridnum,gridnum*4+1)
    num_particles = gridnum*gridnum
    diag = runSolverTest(beam_parameters,mass, num_particles,current_in,gourmet,SOLVER_ORBIT,color='g')
 
    # Gauss solver
    num_particles = gridnum*gridnum
    diag = runSolverTest(beam_parameters,mass, num_particles, current_in,gourmet,SOLVER_GAUSS)
    
    # 3d solver
    num_particles = gridnum*gridnum*gridnum
    #diag = runSolverTest(beam_parameters,mass,num_particles,current_in,gourmet, SOLVER_3D,color='b', periodic=True)
    #diag = runSolverTest(beam_parameters,mass,num_particles,current_in,gourmet, SOLVER_3D,color='k', periodic=False)
    
    # no space charge
    num_particles = gridnum
    diag = runSolverTest(beam_parameters,mass,num_particles,0.,gourmet, NO_SOLVER ,color='r')
    
    pylab.figure(1)
    title = "N=%d, init w = %e, E=%s, I=%s" % (gridnum,xwidth_initial,str(kinetic_energy),str(current_in))
    pylab.title(title)
    
    if do_plot == 1:
	pylab.show()
    #else:
	#pylab.figure(1)
	#filename = 'pic'+str(time.time())+'.png' 
	#pylab.savefig(filename)
    
