import os.path
import sys
import synergia
import Numeric
import numpy
from math import sqrt, sin, acos, pi
import time
import impact

import s2_fish
from mpi4py import MPI

if ( __name__ == '__main__'):

   
    myopts = synergia.Options("PS2")
    myopts.add("latticefile","PS2v.lat","",str)
    myopts.add("maporder",1,"map order",int)
    myopts.add("turns",10,"number of turns",int)
#    myopts.add("latticefile","fobodobo_s.lat","",str)
#    myopts.add("latticefile","fodo.lat","",str)
    myopts.add("xoffset",0.,"transverse offset in x",float)
    myopts.add("yoffset",0.,"transverse offset in y",float)
    myopts.add("emitx",3.24e-06,"X emittance",float)
    myopts.add("emity",1.73e-06,"Y emittance",float)
    myopts.add("sige",1e-3,"(sigma E) over E",float)
    myopts.add("Ekin",4.0,"",float)
    myopts.add("bunchnp",7.0e+11,"number of particles per bunch",float)
    myopts.add("tgridnum",16,"transverse grid cells",int)
    myopts.add("lgridnum",64,"",int)
    myopts.add("bunches",1,"",int)
    myopts.add("partpercell",1,"",float)
    myopts.add("space_charge",1,"",int)
    myopts.add("kicks",40,"kicksper line",int)
    myopts.add("numtrack",10,"number of particles to track",int)
    
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      [myopts.get("latticefile")])
    
    t0 = time.time()
    model="NMCRING"

    
    charge = 1.0  # electron charge in C
    initial_phase = 0.0
#    scaling_frequency = 47713451.5923694
#    scaling_frequency =40.0e6
    scaling_frequency =1.0e9
    
    kinetic_energy= myopts.get("Ekin")
    mass = synergia.PH_NORM_mp
    energy=kinetic_energy+mass
    emitx=myopts.get("emitx")
    emity=myopts.get("emity")
    xoffset = myopts.get("xoffset")
    yoffset = myopts.get("yoffset")
    sige= myopts.get("sige")
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    part_per_cell=myopts.get("partpercell")
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell )
    kicks_per_line = myopts.get("kicks")
    space_charge=myopts.get("space_charge")
    
    
    ee = synergia.Error_eater()
    ee.start()
    
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile"))
        ,model,kinetic_energy,
                        scaling_frequency,myopts.get("maporder"), delay_complete=True)
    
    gourmet.insert_space_charge_markers(kicks_per_line)
    gourmet.complete_setup()
   
    
    
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "line length= ",gourmet.orbit_length()
        #gourmet.print_elements()
   

  
    
    
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    if MPI.COMM_WORLD.Get_rank() ==0:
        print " lattice beta_x, alpha_x, beta_y, alpha_y = ", beta_x, alpha_x, beta_y, alpha_y
    
    
   
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                        initial_phase, scaling_frequency,
                                         transverse=1)	
    beta=beam_parameters.get_beta()
    gamma=beam_parameters.get_gamma() 
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "beta=",beta	
        print "gamma=",gamma			 
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV		
     
    (xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emitx,alpha_x,beta_x)
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
    
    (ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emity,alpha_y,beta_y)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry,offset=yoffset)
       
    bunch_len= 1e-9
    lam_z=sige*gamma*beam_parameters.mass_GeV
    sigma_z_meters = beta*synergia.physics_constants.PH_MKS_c*bunch_len
    if MPI.COMM_WORLD.Get_rank() ==0:
        print " xwidth=",xwidth
	print " ywidth=",ywidth
	print "sigma_z_meters =",sigma_z_meters
	print "sigma z in the code=",beta*synergia.physics_constants.PH_MKS_c/scaling_frequency
	print " "
	print " xpwidth=",xpwidth
	print " ypwidth=",ypwidth   
	print " "
	print " xoffset=",xoffset
	print " yoffset=",yoffset
        
    beam_parameters.z_params(sigma = sigma_z_meters, lam = sige* pz/beta**2)
    # beam_parameters.z_params(sigma = sigma_z_meters, lam = lam_z)
	

    bunch_sp=25.e-9	
    bunch_spacing=beta*synergia.physics_constants.PH_MKS_c*bunch_sp
    
    
    current = myopts.get("bunchnp")* \
        synergia.physics_constants.PH_MKS_e/ \
        (bunch_spacing/(beta*synergia.physics_constants.PH_MKS_c))
    
    if not space_charge:
        current=0.
       
    if MPI.COMM_WORLD.Get_rank() ==0:
       print "Kinetic energy=", kinetic_energy    
       print "bunch_spacing=",bunch_spacing
       print "current =",current
       print "space_charge =",space_charge
       print "num_particles =",num_particles
    
    
    sys.stdout.flush()
     
     
    #BC_choice="trans finite, long periodic round"
    #griddim = (tgridnum+1,tgridnum+1,lgridnum+1)  
    BC_choice="trans open, long periodic"
    griddim = (tgridnum,tgridnum,lgridnum+1)    
    
    num_particles = impact.adjust_particles(
        griddim[0]*griddim[1]*griddim[2] * part_per_cell,MPI.COMM_WORLD.Get_size())
   
	
    
    
    pgrid = impact.Processor_grid(1)
    cgrid = impact.Computational_grid(griddim[0],griddim[1],griddim[2],
                                                  BC_choice)
    piperad =100.
    
    field = impact.Field(beam_parameters, pgrid, cgrid, piperad)
    bunch = impact.Bunch(current, beam_parameters, num_particles, pgrid)
    bunch.generate_particles()          
    bunch.write_particles("beginimpact")
    diag = synergia.Diagnostics_impact(gourmet.get_initial_u())
    
     
    if myopts.get("numtrack") > 0:
        tracker = synergia.Tracker("/tmp",(myopts.get("numtrack"),num_particles))
        tracker.add(bunch,0.0)
    else:
        tracker = None
    track_period_steps=20
    
    
    
    
    log = open("log","w") 
    
    diag1 = synergia.Diagnostics_impact(gourmet.get_initial_u())
    diag1.add(0.0,bunch)
    diag1.write_hdf5("begin_impact")
    
 
    
    
         
    if MPI.COMM_WORLD.Get_rank() ==0:
       output = "start propagation"
       print output
       log.write("%s\n" % output)
       log.flush()
    
    
    s = 0.0  
    for turn in range(1,myopts.get("turns")+1):
       t1 = time.time()
       s = synergia.propagate(s,gourmet,bunch,diag,diag1,griddim,use_impact=True,
        pgrid=pgrid,field=field,cgrid=cgrid,tracker=tracker,track_period_steps=track_period_steps)
       if MPI.COMM_WORLD.Get_rank() ==0:
          output = "turn %d time = %g"%(turn,time.time() - t1)
          print output
          log.write("%s\n" % output)
          log.flush()
    if MPI.COMM_WORLD.Get_rank() == 0:
#           diag.write("mi_impact")
       diag.write_hdf5("mi_impact")
       diag1.write_hdf5("bk_impact")  
    if tracker:
      tracker.close()
      tracker.show_statistics()   
	  	  
    log.close()
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "elapsed time =",time.time() - t0
    
