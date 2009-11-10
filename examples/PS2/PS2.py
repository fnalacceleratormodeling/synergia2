import os.path
import sys
import synergia
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
#    myopts.add("xoffset",3.e-7,"transverse offset in x",float)
#    myopts.add("yoffset",3.e-7,"transverse offset in y",float)
    myopts.add("xoffset",0.000001,"transverse offset in x",float)
    myopts.add("yoffset",0.000001,"transverse offset in y",float)
    #myopts.add("emitx",3e-06,"X emittance",float)
    #myopts.add("emity",3e-06,"Y emittance",float)
    #myopts.add("emitz",0.098,"Z emittance",float)
    #myopts.add("sige",10.e-3,"(sigma E) over E",float)
    #myopts.add("sige",1.e-3,"(sigma E) over E",float)
    myopts.add("xrms",4.51e-3," xrms", float)
    myopts.add("yrms",2.81e-3," yrms", float)
    myopts.add("trms_rad",1.11," trms in radians, bucket size=2*pi", float)
    myopts.add("Ekin",4.0,"",float)
    myopts.add("bunchnp",4.2e+11,"number of particles per bunch",float)
    myopts.add("tgridnum",32,"transverse grid cells",int)
    myopts.add("lgridnum",32,"",int)
   # myopts.add("bunches",1,"",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("outputperiod",50,"save beam each <saveperiod> turns",int)
    myopts.add("partpercell",2,"",float)
    myopts.add("space_charge",0,"",int)
    myopts.add("kicks",60,"kicksper line",int)
    myopts.add("numtrack",0,"number of particles to track",int)
    myopts.add("solver","3dc","solver type for spch",str)
    myopts.add("exactrf",1,"use exact propagator for RF cavities",int)
    myopts.add_suboptions(synergia.opts)
    myopts.parse_argv(sys.argv)
    job_mgr = synergia.Job_manager(sys.argv,myopts,
                                      [myopts.get("latticefile")])
    
    
    t0 = time.time()
    ee = synergia.Error_eater()
    ee.start()
# making gourmet    

    model="NMCRING"
    kinetic_energy= myopts.get("Ekin")
    mass = synergia.PH_NORM_mp
    energy=kinetic_energy+mass
    #scaling_frequency=39349217.5333
    scaling_frequency=40e6
    
    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile"))
        ,model,kinetic_energy,
                        scaling_frequency,myopts.get("maporder"), delay_complete=True,
                        exact_rf=myopts.get("exactrf"))
    
    
    kicks_per_line = myopts.get("kicks")
    gourmet.insert_space_charge_markers(kicks_per_line)
    
    
    gamma = energy/mass
    beta = sqrt(1.0-1.0/gamma**2)
    w_cav=beta*synergia.physics_constants.PH_MKS_c/gourmet.orbit_length()
    print "w_cav=",w_cav," w_cav*180=",w_cav*180
 
    for element in gourmet.beamline:    
        if element.Type() == 'rfcavity':
            element.setFrequency(scaling_frequency)
            element.setPhi(0.)
           # strength = 10.e-9*1.5e6
            #element.setStrength(strength)
            
    gourmet.complete_setup()
    
 
    
    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    (tune_x, tune_y, tune_z)           = synergia.matching.get_tunes(gourmet)
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
        print "(tune_x, tune_y, tune_z) = %g, %g,  %g" % (tune_x, tune_y, tune_z)
    
    
 # defining beam_parameters
    charge = 1.0
    initial_phase = 0.0
    
    beam_parameters = synergia.Beam_parameters(mass, charge, kinetic_energy,
                                         initial_phase, scaling_frequency_Hz=scaling_frequency,
                                         transverse=0, adjust_zlength_to_freq=1)
    
    beta=beam_parameters.get_beta()
    gamma=beam_parameters.get_gamma() 
    pz = gamma * beta * beam_parameters.mass_GeV
    
   
    
    full_map = gourmet.get_single_linear_map()   
    xwidth= myopts.get("xrms")
    ywidth= myopts.get("yrms")
    zwidth= (myopts.get("trms_rad")/(2*pi))*beam_parameters.get_z_length()
   

    rms_index=[0,2,4]
    C=synergia.rms_match_3d(full_map[0:6,0:6],beam_parameters,xwidth,ywidth,zwidth,rms_index,print_emittances=True)
    beam_parameters.offset_x_m = myopts.get("xoffset")
    beam_parameters.offset_y_m = myopts.get("yoffset")
    
    
    if MPI.COMM_WORLD.Get_rank() ==0:
        print " **********************************************************************" 
        print "beam information:"
        print
        print "transverse=",beam_parameters.transverse
        if beam_parameters.adjust_zlength_to_freq:
            print "adjusted length =",beam_parameters.z_length
        else:
            print "bunch length (not adjusted)= ", beam_parameters.z_length
           
    sys.stdout.flush()
    
    
    
    
    
    
### creating the bunch 
    bunchnp=myopts.get("bunchnp") 
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    part_per_cell = myopts.get("partpercell")
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * myopts.get("partpercell"))
    diag=synergia.Diagnostics(gourmet.get_initial_u(),short=True)
    #bunch = s2_fish.Macro_bunch.gaussian(bunchnp,num_particles,beam_parameters,diagnostics=diag,periodic=True)
    bunch=s2_fish.Macro_bunch.gaussian_covariance(bunchnp,num_particles,beam_parameters,C,diagnostics=diag,periodic=True)    
    bunch.write_particles("begin_particles")  
  
    
    
    
    #pipe_radius=0.08
    #if space_charge:
        #solver="s2_fish_cylindrical"
       ## griddim = (16,16,33)
       ## solver="s2_fish_3d"

        #sp_ch=s2_fish.SpaceCharge(solver,griddim,radius_cylindrical=pipe_radius,periodic=True) 
        #if MPI.COMM_WORLD.Get_rank() ==0:  
            #print " sp_ch grid=",sp_ch.get_grid()
            #print " sp_ch solver=",sp_ch.get_solver()
            #print " sp_ch pipe radius=",sp_ch.get_radius_cylindrical()
    #else:
       #sp_ch=None  
    
    space_charge=myopts.get("space_charge")
    if space_charge:
        solver=myopts.get("solver")
        if (solver=="3DC") or (solver=="3dc"):
            griddimsp = (2*tgridnum,tgridnum,lgridnum)
            solversp="s2_fish_cylindrical"
            spch_pipe_radius = 0.08
            sp_ch=s2_fish.SpaceCharge(solversp,griddimsp,radius_cylindrical=spch_pipe_radius,periodic=True) 
            if MPI.COMM_WORLD.Get_rank() ==0: 
                print " sp_ch grid=",sp_ch.get_grid()
                print " sp_ch solver=",sp_ch.get_solver()
                print " sp_ch pipe radius=",sp_ch.get_radius_cylindrical()
        elif (solver=="3D") or (solver=="3d"): 
                griddimsp = (tgridnum,tgridnum,4*lgridnum+1)        
                solversp="s2_fish_3d"
                sp_ch=s2_fish.SpaceCharge(solversp,grid=griddimsp,periodic=True)
                if MPI.COMM_WORLD.Get_rank() ==0: 
                    print " sp_ch grid=",sp_ch.get_grid()
                    print " sp_ch solver=",sp_ch.get_solver()
        else:
             raise RuntimeError,  " Choose either a 3d or a 3dc solver "           
    else:
       sp_ch=None  
    
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "macroparticles =",num_particles
        print "bunchnp = ", myopts.get("bunchnp")
        print "maporder = ", myopts.get("maporder")
        print "kicks/line = ", myopts.get("kicks")
        print "offsets x,y,z: ", myopts.get("xoffset"), \
            myopts.get("yoffset"), myopts.get("zoffset")
        print "emittance: ", myopts.get("emittance")
        #print "bunch length = ", myopts.get("bunchlen")
      #  print "dpop width = ", myopts.get("dpop")
        print "using lattice file: ", myopts.get("latticefile")
        print "line length= ",gourmet.orbit_length()
        print "beta=",beta  
        print "gamma=",gamma   
        print "beta*gamma=",beta*gamma   
          
   
    sys.stdout.flush()
  
    
  
 

    
    if myopts.get("numtrack") > 0:
        tracker = synergia.Tracker("/tmp",(myopts.get("numtrack"),num_particles))
        tracker.add(bunch,0.0)
    else:
        tracker = None
    track_period_steps=20
    
    
    
    log = open("log","w") 
    if MPI.COMM_WORLD.Get_rank() ==0:
       output = "start propagation"
       print output
       log.write("%s\n" % output)
       log.flush()
    
#start propagation  
    outputperiod=myopts.get("outputperiod")
    s=0.0
    for turn in range(1,myopts.get("turns")+1):
       t1 = time.time()
       if turn % myopts.get("saveperiod") == 0:
            bunch.write_particles("turn_%03d.h5" %(turn-1))
       
       s = synergia.propagate(s,gourmet, bunch,  space_charge=sp_ch)
       if turn % outputperiod==0:   
            bunch.diagnostics.write_hdf5("tmp_output-%02d" %(turn/outputperiod)) 
      # bunch.write_particles("turn_%03d.h5" % turn)
     
       if MPI.COMM_WORLD.Get_rank() ==0:
          output = "turn %d time = %g"%(turn,time.time() - t1)
          print output
          log.write("%s\n" % output)
          log.flush()
    

    bunch.diagnostics.write_hdf5("ps2_output")
    bunch.write_particles("end_particles")       
    

    if tracker:
       tracker.close()
       tracker.show_statistics()    
    
    log.close()
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "elapsed time =",time.time() - t0
    
