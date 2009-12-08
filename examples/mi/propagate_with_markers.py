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

    myopts = synergia.Options("mi")
    myopts.add("transverse",0,"longitudinally uniform beam",int)
    myopts.add("maporder",2,"map order",int)
    myopts.add("emittance",6.032e-06,"emittance",float)
    # monochromatic beam for now
    myopts.add("dpop",1.0e-20,"(delta p)/p RMS width",float)
    myopts.add("bunchlen", 0.5, "RMS bunchs length (z width) [m]", float)
    myopts.add("dpopoffset", 0.0, "offset in dpop", float)
    myopts.add("kicks",208,"kicks per line",int)
    myopts.add("turns",10,"number of turns",int)
    myopts.add("latticefile","mi20-egs.lat","",str)
    myopts.add("tgridnum",16,"transverse grid cells",int)
    myopts.add("lgridnum",64,"",int)
    myopts.add("xoffset",1.0e-6,"transverse offset in x",float)
    myopts.add("yoffset",-2.0e-6,"transverse offset in y",float)
    myopts.add("xpoffset", 0, "offset in x-prime", float)
    myopts.add("ypoffset", 0, "offset in y-prime", float)
#    myopts.add("zoffset",0,"offset in z", float)
#    myopts.add("xoffset",4.26e-4,"transverse offset in x",float)
#    myopts.add("yoffset",1.86e-4,"transverse offset in y",float)
    myopts.add("zoffset",0,"offset in z", float)
    myopts.add("space_charge",1,"apply space charge kicks",int)
    myopts.add("impedance",0,"apply resistive wall impedance kicks",int)
    myopts.add("pipe_symmetry","circular","",str) 
  #  myopts.add("pipe_symmetry","x_parallel_plates","",str) 
    myopts.add("energy",8.9,"",float)
    myopts.add("partpercell",1,"",float)
    myopts.add("bunches",1,"",int)
    myopts.add("bunchnp",1.0e11,"number of particles per bunch",float)
    myopts.add("partdump", 0, "dump particles each turn", int)
    
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      [myopts.get("latticefile")])

    t0 = time.time()
   
    ee = synergia.Error_eater()
    ee.start()
    
    
### making gourmet 
    energy = myopts.get("energy")
    mass = synergia.PH_NORM_mp
    kinetic_energy = energy-mass
    # this scaling frequency sets the length unit scaling to unity
    scaling_frequency = synergia.PH_MKS_c/(2.0*math.pi)
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile"))
        ,"ring_p_q605",kinetic_energy,
                        scaling_frequency, myopts.get("maporder"),delay_complete=True)
	
    # introduce kicks
    kicks_per_line = myopts.get("kicks")		
    gourmet.insert_space_charge_markers(kicks_per_line) 		

    newstdout = open("mi20-egs-with-markers.txt","w")
    sys.stdout.flush()
    oldstdout = sys.stdout
    sys.stdout = newstdout
    print "hi there"
    gourmet.print_elements()
    sys.stdout.flush()
    sys.stdout = oldstdout

    # set rf freq  
    line_length = gourmet.orbit_length()
    gamma = energy/mass
    beta = math.sqrt(1.0-1.0/gamma**2)

    # revolution cavity frequency
    w_cav=beta*synergia.physics_constants.PH_MKS_c/line_length

    if MPI.COMM_WORLD.Get_rank() == 0:
	    print "line length: ", line_length
	    print "gamma initial: ",gamma
	    print "1-beta initial: ", 1-beta
	    print "revolution frequency (w_cav): ", w_cav
	    sys.stdout.flush()
	    
    	
    
    # walk through the beamline and activate RF cavities
    activate_rf_cavities = False
    
    if activate_rf_cavities:

        for element in gourmet.beamline:	
            if element.Type() == 'rfcavity':
                #	   element.setHarmonicNumber(4)	
                #	   h=element.getHarmonicNumber()
                # print"h=",h	  
                #element.setFrequency(59955852.5381452)
                element.setFrequency(32*w_cav)
                # print "cavitiy w=",element.getRadialFrequency()/(2.*math.pi)
                element.setPhi(math.pi)


    gourmet.complete_setup()

    sys.stdout.flush()
    oldstdout = sys.stdout
    newstdout = open("beamline_trace.txt", "w")
    sys.stdout = newstdout

    testpart = gourmet.get_initial_particle()
    print "Start of beamline trace, particle: ", gourmet.printpart(testpart)
    print "%5s %10s %5s %13s %s: |x npx y npy cdt ndp|" % ("index", 's_begin', 'length', 'type','name')

    s = 0.0
    idx = 0
    for elem in gourmet.beamline:
        lenelem = elem.OrbitLength(gourmet.get_initial_particle())
        elem.propagate(testpart)
        print "%5s %10s %5s %13s %s: " % \
              (idx, s, lenelem, elem.Type(), elem.Name()),
        idx += 1
        s += lenelem
        gourmet.printpart(testpart)
        if (abs(testpart.get_x()) > 1.0e-6) or \
           (abs(testpart.get_npx()) > 1.0e-6) or \
           (abs(testpart.get_y()) > 1.0e-6) or \
           (abs(testpart.get_npy()) > 1.0e-6) or \
           (abs(testpart.get_cdt()) > 1.0e-6) or \
           (abs(testpart.get_ndp()) > 1.0e-6):
            print "!!!  BAD! BAD! BAD!  particle deviating from 0!!!"
        sys.stdout.flush()

    print "End of line, s = %10s" % s
    sys.stdout.flush()
    sys.stdout = oldstdout
    sys.exit(20)

    print "checking gourmet"
    gourmet.check()
    sys.stdout.flush()
  
   
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    (tune_x, tune_y, tune_z)           = synergia.matching.get_tunes(gourmet)
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
        print "(tune_x, tune_y, tune_z) = %g, %g,  %g" % (tune_x, tune_y, 1.-tune_z)
	print "initial_u: ", gourmet.get_initial_u()
	sys.stdout.flush()

    
# defining beam_parameters
    charge = 1.0
    initial_phase = 0.0
    
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency_Hz=scaling_frequency,
                                         transverse=0, adjust_zlength_to_freq=1)
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
    
    
    zoffset = myopts.get("zoffset")
    zwidth=myopts.get("bunchlen")
    zpwidth=myopts.get("dpop")
    beam_parameters.z_params(sigma = zwidth, 
                             lam = zpwidth* pz, z_length=0.33, offset=zoffset,
                             offset_p = myopts.get("dpopoffset")*pz)
  
   
    sys.stdout.flush()

 

### creating the bunch
    bunchnp=myopts.get("bunchnp") 
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    part_per_cell = myopts.get("partpercell")
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell)
    diag=synergia.Diagnostics(gourmet.get_initial_u(),short=True)
    bunch = s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,periodic=True,diagnostics=diag)
    #bunch= s2_fish.Macro_bunch.from_bunch(bunch1)
    #bunch = s2_fish.Macro_bunch.test(int(part_per_cell))
   # bunch = s2_fish.Macro_bunch.test_am(bunchnp,part_per_cell,griddim,beam_parameters)
   # bunch = s2_fish.Macro_bunch.sphere(num_particles, 0.001)
   # bunch = s2_fish.Macro_bunch.cylinder(num_particles, 0.001,0.01)
   # covar=beam_parameters.get_covariances()
   # bunch = s2_fish.Macro_bunch.gaussian_covariance(bunchnp,num_particles,beam_parameters,covar,periodic=True)
   # print "bunch first long size=",bunch.get_longitudinal_period_size()
   
    log = open("log","w")
    if MPI.COMM_WORLD.Get_rank() ==0:
            output = "start propagation"
            print output
            log.write("%s\n" % output)
            log.flush()

 
 
    pipe_radius=0.025
    
    
    space_charge=myopts.get("space_charge")
    if space_charge:
        #griddim = (16,16,32)
        #solver="s2_fish_cylindrical"
        griddim = (tgridnum,tgridnum,lgridnum)
        solver="s2_fish_3d"

        sp_ch=s2_fish.SpaceCharge(solver,griddim,radius_cylindrical=pipe_radius,periodic=True)	
        print " sp_ch grid=",sp_ch.get_grid()
        print " sp_ch solver=",sp_ch.get_solver()
        print " sp_ch pipe radius=",sp_ch.get_radius_cylindrical()
    else:
	   sp_ch=None    	




 
    impedance=myopts.get("impedance")
    if impedance:
        pipe_conduct= 1.4e6 # [ohm^-1 m^-1] (stainless steel)
        wall_thickness=0.0114        
        pipe_symmetry=myopts.get("pipe_symmetry")
        rw_impedance=s2_fish.Impedance(pipe_radius, pipe_conduct,wall_thickness, line_length,lgridnum, pipe_symmetry)
                #pipe_symmetry="x_parallel_plates")
        print "IMPEDANCE PIPE radius=", rw_impedance.get_pipe_radius()
        print "IMPEDANCE PIPE wall_thickness=",rw_impedance.get_wall_thickness()
        print "IMPEDANCE PIPE symmetry=",rw_impedance.get_pipe_symmetry()
    else:
        rw_impedance=None      
    
    
    
    
    if MPI.COMM_WORLD.Get_rank() ==0:
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
        print "sp_ch = ", sp_ch
    
    # write out initial particle distribution
    bunch.write_particles("begin-00")
    if myopts.get("partdump") != 0:
        bunch.write_particles("turn-%04d" % 0)
    
### begin propagation    
    s = 0.0
    for turn in range(1,myopts.get("turns")+1):
        t1 = time.time()
        s = synergia.propagate(s,gourmet, bunch, space_charge=sp_ch,impedance=rw_impedance)
	    	
        if MPI.COMM_WORLD.Get_rank() ==0:
            output = "turn %d time = %g"%(turn,time.time() - t1)
            print output
            log.write("%s\n" % output)
            log.flush()

        # dump particles after this turn if requested
        if myopts.get("partdump") != 0:
	    #sys.stdout.write("dumping particles for turn %04d\n" %  turn)
	    #sys.stdout.flush()
	    #log.write("dumping particles for turn %04d\n" %  turn)
	    #log.flush()
            bunch.write_particles("turn-%04d" % turn)

    sys.stdout.write("Writing saved diagnostics\n")
    sys.stdout.flush()

    log.write("Writing saved diagnostics\n")
    log.flush()

    bunch.diagnostics.write_hdf5("mi-00")
    
    sys.stdout.write("Writing end particles\n")
    sys.stdout.flush()

    log.write("Writing end particles")
    log.flush()

    # write final particle distribution
    bunch.write_particles("end-00")

    log.close()
    
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "elapsed time =",time.time() - t0
 
