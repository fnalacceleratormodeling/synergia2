#!/usr/bin/env bwpython

import numpy
import time
import math
import os
import sys

import synergia
import s2_fish

from mpi4py import MPI

write_memory_diagnostics_by_rank = False  # write memory usage each rank
write_memory_diagnostics_to_output = True # write memory usage to stdout
write_memory_diagnostics_each = 1000 # write memory usage out for every n turns

if ( __name__ == '__main__'):
    t0 = time.time()

    myopts = synergia.Options("circular")
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",1,"map order",int)
    myopts.add("emittance",5.89533703303356e-07,"emittance",float)
    # longitudinal beta is 143.6
    myopts.add("dpop",3.482e-4,"(delta p)/p RMS width",float)
    myopts.add("bunchlen", 0.05, "RMS bunchs length (z width) [m]", float)
    myopts.add("dpopoffset", 0.0, "offset in dpop", float)
    myopts.add("kicks",32,"kicks per line",int)
    myopts.add("turns",10,"number of turns",int)
    myopts.add("latticefile","foborodobo_s.lat","",str)
    myopts.add("tgridnum",16,"transverse grid cells",int)
    myopts.add("lgridnum",64,"",int)
    myopts.add("xoffset",0.0004,"transverse offset in x",float)
    myopts.add("yoffset",-0.0002,"transverse offset in y",float)
    myopts.add("xpoffset", 0, "offset in x-prime", float)
    myopts.add("ypoffset", 0, "offset in y-prime", float)
#    myopts.add("zoffset",0,"offset in z", float)
#    myopts.add("xoffset",4.26e-4,"transverse offset in x",float)
#    myopts.add("yoffset",1.86e-4,"transverse offset in y",float)
    myopts.add("zoffset",0,"offset in z", float)
    myopts.add("space_charge",0,"",int)
    myopts.add("impedance",0,"",int)
    myopts.add("energy",100.004401675138,"",float)
    myopts.add("partpercell",1,"",float)
    myopts.add("bunches",1,"",int)
    myopts.add("bunchnp",1.0e11,"number of particles per bunch",float)
    
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
    #scaling_frequency = 47713451.5923694
    #scaling_frequency = 1000000
    scaling_frequency = 1
    pipexradius = 0.03
    pipeyradius = 0.03
#    pipexradius = 0.123
#    pipeyradius = 0.0508
    part_per_cell = myopts.get("partpercell")
    kicks_per_line = myopts.get("kicks")
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell)
    solver = "3d"
    xoffset = myopts.get("xoffset")
    xpoffset = myopts.get("xpoffset")
    yoffset = myopts.get("yoffset")
    ypoffset = myopts.get("ypoffset")
    zoffset = myopts.get("zoffset")
    
    pipe_conduct= 1.4e6 # [ohm^-1 m^-1] (stainless steel)
    
    impedance=myopts.get("impedance")
    space_charge=myopts.get("space_charge")

    if MPI.COMM_WORLD.Get_rank() ==0:
        print "space_charge =",space_charge
        print "impedance =", impedance
        print "macroparticles =",num_particles
        print "bunchnp = ", myopts.get("bunchnp")
        print "maporder = ", myopts.get("maporder")
        print "kicks/line = ", myopts.get("kicks")
        print "offsets x,y,z: ", myopts.get("xoffset"), \
              myopts.get("yoffset"), myopts.get("zoffset")
        print "emittance: ", myopts.get("emittance")
        print "bunch length = ", myopts.get("bunchlen")
        print "dpop width = ", myopts.get("dpop")
        print "using lattice file: ", myopts.get("latticefile")
        print "grid = ", myopts.get("tgridnum"),"^2 x ",\
              myopts.get("lgridnum")

    ee = synergia.Error_eater()
    ee.start()
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile"))
        ,"model",kinetic_energy,
                        scaling_frequency, myopts.get("maporder"),delay_complete=True)
			
    gourmet.insert_space_charge_markers(kicks_per_line) 		


### set rf freq  
    line_length = gourmet.orbit_length()
    gamma = energy/mass
    beta = math.sqrt(1.0-1.0/gamma**2)



    w_cav=beta*synergia.physics_constants.PH_MKS_c/line_length
    print "w_cav=",w_cav
    	
    
    for element in gourmet.beamline:	
	if element.Type() == 'rfcavity':
#	   element.setHarmonicNumber(4)	
#	   h=element.getHarmonicNumber()
	   #print"h=",h	  
	   #element.setFrequency(59955852.5381452)
	   element.setFrequency(32*w_cav)
	  # print "cavitiy w=",element.getRadialFrequency()/(2.*math.pi)
	   element.setPhi(math.pi)



    gourmet.complete_setup()






   

    # try without commissioning
    #gourmet.needs_commission = True
    #gourmet.is_commissioned = False
   # gourmet._commission()

    #gourmet.check(1)
    #jet_particle = gourmet.get_initial_jet_particle()
    #particle = gourmet.get_initial_particle()
    
  
    #gourmet.printpart(particle)
    #for element in gourmet.beamline:
	#energy = jet_particle.ReferenceEnergy()
	#print element.Name(),element.Type()
	#jet_particle = gourmet.get_jet_particle(energy)
	#print "xxxxxxxxxxx"
	#print element.Name(),element.Type()
	#gourmet.printpart(particle)
	
        #if element.Type() == 'rfcavity':
	   #print "before"
	   #print "particle state before=",particle.State()
           #element.propagate(jet_particle)
           #element.propagate(particle)
	  
	   
           #map = gourmet._convert_linear_maps([jet_particle.State().jacobian()])[0]
	   #print numpy.array2string(map,precision=1)
	   #print "after"
	   #print "particle state after=",particle.State()
	   #gourmet.printpart(particle)
   
   
   
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    (tune_x, tune_y)                    = synergia.matching.get_tunes(gourmet)
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
        print "(tune_x, tune_y) = %g, %g" % (tune_x, tune_y)
	
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency,
                                         transverse=0)
    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 

    emittance = myopts.get("emittance")
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV

    (xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emittance,alpha_x,beta_x)
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,
                             r = rx,offset=xoffset, offset_p = xpoffset * pz)
    
    (ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emittance,alpha_y,beta_y)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,
                             r = ry,offset=yoffset, offset_p = ypoffset * pz )
    
    beam_parameters.z_params(sigma = myopts.get("bunchlen"),
                             lam = myopts.get("dpop")* pz, offset=zoffset,
                             offset_p = myopts.get("dpopoffset")*pz)

    sys.stdout.flush()

    s = 0.0
    bunch_spacing = line_length/myopts.get("bunches")
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "line_length =",line_length
        print "bunch_spacing =",bunch_spacing
    tau = 0.5*line_length/kicks_per_line
    kick_time = 0.0
    beta = beam_parameters.get_beta()
    # current is (bunch charge)/(1 period of scaling frequency)
    if not (space_charge or impedance) :
	current =0.
    else:
         current = myopts.get("bunchnp") * synergia.physics_constants.PH_MKS_e * \
              scaling_frequency
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "current =",current
         
        

   # numbunches = myopts.get("bunches")
   # bunches = []
   # diags = []
    #for bunchnum in range(0,numbunches):
        #bunches.append(s2_fish.Macro_bunch(mass,1))
        #bunches[bunchnum].init_gaussian(num_particles,current,beam_parameters)
        #bunches[bunchnum].write_particles("begin-%02d"%bunchnum)
        #diags.append(synergia.Diagnostics(gourmet.get_initial_u(),short=True))


    bunch = s2_fish.Macro_bunch(mass,1)
    bunch.init_gaussian(num_particles,current,beam_parameters)    
    bunch.write_particles("begin")
    diag = synergia.Diagnostics(gourmet.get_initial_u(),save_period=0)   

    log = open("log","w")
    if MPI.COMM_WORLD.Get_rank() ==0:
            output = "start propagation"
            print output
            log.write("%s\n" % output)
            log.flush()


    if impedance:
        rw_impedance=s2_fish.Impedance(pipexradius, pipe_conduct, line_length,lgridnum, pipe_symmetry="circular")
	#pipe_symmetry="x parallel plates")
        print "IMPEDANCE PIPE radius=", rw_impedance.get_pipe_radius()
    else:
	rw_impedance=None      
    
    
    
    for turn in range(1,myopts.get("turns")+1):
        t1 = time.time()

        s = synergia.propagate(s,gourmet,
            bunch,diag,griddim,use_s2_fish=True,periodic=True,
            space_charge=space_charge, rw_impedance=rw_impedance)
	    	
        if MPI.COMM_WORLD.Get_rank() ==0:
            output = "turn %d time = %g"%(turn,time.time() - t1)
            print output
            log.write("%s\n" % output)
            log.flush()

   
    if MPI.COMM_WORLD.Get_rank() == 0:
        diag.write_hdf5("mi-00")
    
    bunch.write_particles("end-00")

    if MPI.COMM_WORLD.Get_rank() ==0:
        print "elapsed time =",time.time() - t0
 
 
