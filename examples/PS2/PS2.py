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

def ha_match(map,beam_parameters,emitx,emity,dpop):
    numpy_map = map
    evals,evect_matrix = numpy.linalg.eig(numpy_map)
    evects = []
    for i in range(0,6):
        evects.append(evect_matrix[:,i])
    E = range(0,3)
    remaining = range(5,-1,-1)
    for i in range(0,3):
        # find complex conjugate among remaining eigenvectors
        first = remaining.pop()
        best = 1.0e30
        conj = -1
        for item in remaining:
            sum = evects[first]+evects[item]
            if abs(numpy.max(sum.imag)) < best:
                best = abs(numpy.max(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError,"failed to find a conjugate pair in ha_match"
        remaining.remove(conj)
        #~ print first,conj,best
        tmp=numpy.outer(evects[first],
            numpy.conjugate(evects[first]))
        tmp+=numpy.outer(evects[conj],
            numpy.conjugate(evects[conj]))
        E[i]=tmp.real

    for srt in range(1,3):
        if abs(numpy.linalg.det(E[srt][0:2,0:2])) > abs(numpy.linalg.det(E[0][0:2,0:2])):
            tmp = E[0]
            E[0] = E[srt]
            E[srt] = tmp
    if abs(numpy.linalg.det(E[2][2:4,2:4])) > abs(numpy.linalg.det(E[1][2:4,2:4])):
        tmp = E[1]
        E[1] = E[2]
        E[2] = tmp
    Cxy, Cxpyp, Cz, Czp = beam_parameters.get_conversions()
    gamma = beam_parameters.get_gamma()
    C = numpy.zeros([6,6],'d')
    C += E[0]*emitx/(Cxpyp*sqrt(abs(numpy.linalg.det(E[0][0:2,0:2]))))
    C += E[1]*emity/(Cxpyp*sqrt(abs(numpy.linalg.det(E[1][2:4,2:4]))))
    C += E[2]*dpop*dpop/(gamma*gamma*E[2][5,5])
    
    return C

def ha_matche(map,beam_parameters,emitx,emity,emitz):
    numpy_map = map
    evals,evect_matrix = numpy.linalg.eig(numpy_map)
    evects = []
    for i in range(0,6):
        evects.append(evect_matrix[:,i])
    E = range(0,3)
    remaining = range(5,-1,-1)
    for i in range(0,3):
        # find complex conjugate among remaining eigenvectors
        first = remaining.pop()
        best = 1.0e30
        conj = -1
        for item in remaining:
            sum = evects[first]+evects[item]
            if abs(numpy.max(sum.imag)) < best:
                best = abs(numpy.max(sum.imag))
                conj = item
        if conj == -1:
            raise RuntimeError,"failed to find a conjugate pair in ha_match"
        remaining.remove(conj)
        #~ print first,conj,best
        tmp=numpy.outer(evects[first],
            numpy.conjugate(evects[first]))
        tmp+=numpy.outer(evects[conj],
            numpy.conjugate(evects[conj]))
        E[i]=tmp.real

    for srt in range(1,3):
        if abs(numpy.linalg.det(E[srt][0:2,0:2])) > abs(numpy.linalg.det(E[0][0:2,0:2])):
            tmp = E[0]
            E[0] = E[srt]
            E[srt] = tmp
    if abs(numpy.linalg.det(E[2][2:4,2:4])) > abs(numpy.linalg.det(E[1][2:4,2:4])):
        tmp = E[1]
        E[1] = E[2]
        E[2] = tmp
    Cxy, Cxpyp, Cz, Czp = beam_parameters.get_conversions()
    gamma = beam_parameters.get_gamma()
    C = numpy.zeros([6,6],'d')
    C += E[0]*emitx/(Cxpyp*sqrt(abs(numpy.linalg.det(E[0][0:2,0:2]))))
    C += E[1]*emity/(Cxpyp*sqrt(abs(numpy.linalg.det(E[1][2:4,2:4]))))
    C += E[2]* emitz/(Czp* sqrt(abs(numpy.linalg.det(E[2][4:6,4:6]))))    
    #C += E[2]*dpop*dpop/(gamma*gamma*E[2][5,5])
    
    return C

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
    myopts.add("emitx",3e-06,"X emittance",float)
    myopts.add("emity",3e-06,"Y emittance",float)
    myopts.add("emitz",0.098,"Z emittance",float)
    myopts.add("sige",10.e-3,"(sigma E) over E",float)
    #myopts.add("sige",1.e-3,"(sigma E) over E",float)
    myopts.add("Ekin",4.0,"",float)
    myopts.add("bunchnp",4.2e+11,"number of particles per bunch",float)
    myopts.add("tgridnum",32,"transverse grid cells",int)
    myopts.add("lgridnum",32,"",int)
   # myopts.add("bunches",1,"",int)
    myopts.add("saveperiod",10,"save beam each <saveperiod> turns",int)
    myopts.add("outputperiod",100,"save beam each <saveperiod> turns",int)
    myopts.add("partpercell",2,"",float)
    myopts.add("space_charge",0,"",int)
    myopts.add("kicks",40,"kicksper line",int)
    myopts.add("numtrack",0,"number of particles to track",int)
    myopts.add("solver","3dc","solver type for spch",str)
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
                        scaling_frequency,myopts.get("maporder"), delay_complete=True)
    
    
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
    
    #for element in gourmet.beamline:
        #if element.Type() == 'rfcavity':
            #print " rf frecv=", element.getRadialFrequency()
    
    
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
    
    #emitx=myopts.get("emitx")
    #(xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emitx,alpha_x,beta_x)
    #if MPI.COMM_WORLD.Get_rank() == 0:
        #print "xwidth=",xwidth
    #xoffset = myopts.get("xoffset")
    #beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
    
    
    #emity=myopts.get("emity")
    #yoffset = myopts.get("yoffset")  
    #(ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emity,alpha_y,beta_y)
    #if MPI.COMM_WORLD.Get_rank() == 0:
        #print "ywidth=",ywidth
    #beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry,offset=yoffset)
    
    #sige= myopts.get("sige")
    #bunch_len= 1e-9
    #lam_z=sige*energy
    
    #sigma_z_meters = beam_parameters.get_z_length()/6.#beta*synergia.physics_constants.PH_MKS_c*bunch_len
    #print "lam_z=",lam_z
    #beam_parameters.z_params(sigma = sigma_z_meters, lam = lam_z, offset=0.00001)
    
    
    full_map = gourmet.get_single_linear_map()    
    C = ha_match(full_map[0:6,0:6],beam_parameters,myopts.get("emitx"),
                 myopts.get("emity"),myopts.get("sige")/(beta*beta))
    beam_parameters.offset_x_m = myopts.get("xoffset")
    beam_parameters.offset_y_m = myopts.get("yoffset")
   
                 
                 
    
    if MPI.COMM_WORLD.Get_rank() ==0:
        (Cxy, Cxpyp, Cz, Czp)=beam_parameters.get_conversions()
        print "xwidth=",sqrt(C[0,0])/Cxy
        print "ywidth=",sqrt(C[2,2])/Cxy
        print "twidth=",sqrt(C[4,4])/(Cz*synergia.physics_constants.PH_MKS_c*beta)
        print "zmrs=",sqrt(C[4,4])/Cz    
        print "twidth in rad=", sqrt(C[4,4])/Cz*2.*pi/gourmet.orbit_length()
        print "revol time=", gourmet.orbit_length()/(beta*synergia.physics_constants.PH_MKS_c)
        print "bucket time=",1./scaling_frequency
        print "bucket in rad=", (1./scaling_frequency)*beta*synergia.physics_constants.PH_MKS_c*2.*pi/gourmet.orbit_length()
        print "bucket length=",beam_parameters.get_z_length()
        print "lam_z=",sqrt(C[5,5]/Czp)
    
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
    
