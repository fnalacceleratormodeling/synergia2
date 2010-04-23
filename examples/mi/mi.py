#!/usr/bin/env python

import numpy
import time
import math
import os
import sys

import synergia
import s2_fish
import physics_toolkit

from mi_options import myopts

from mpi4py import MPI

DEBUG = False

if ( __name__ == '__main__'):
    t0 = time.time()


    t0 = time.time()

    # get options relating to checkpointing
    do_checkpoint = myopts.get("checkpoint")
    if do_checkpoint:
        checkpoint_freq = myopts.get("checkpoint_freq")
        # the name that will be used to save particle data.  In principle
        # this could be made an option.  The actual name will be
        # ckptsave-bbb-rrrrr.h5 where bbb is the bunch number and rrrrr is the
        # MPI processor rank.
        checkpoint_name = "ckptsave"

    do_resume = myopts.get("resumejob")
    resume_dir = myopts.get("resumedir")
    
    do_dump = myopts.get("dump")

    energy = myopts.get("energy")
    mass = synergia.PH_NORM_mp
    kinetic_energy = energy-mass
    charge = 1.0
    initial_phase = 0.0
    scaling_frequency = myopts.get("scaling_frequency")
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
    if DEBUG:
        print "before Gourmet"
        sys.stdout.flush()

    gourmet = synergia.Gourmet(os.path.join(os.getcwd(),myopts.get("latticefile")),"ring_p_q605",
                               kinetic_energy,
                               scaling_frequency,
                               order=myopts.get("maporder"),
                               delay_complete=True)

    if DEBUG:
        print "after Gourmet"
        sys.stdout.flush()

    # slots at the ends of drifts to match bends
    armando = physics_toolkit.DriftConverter()
    if DEBUG:
        print "after armando"
        sys.stdout.flush()

    gourmet.beamline = armando.convert(gourmet.beamline)
    if DEBUG:
        print "after armando.convert"
        sys.stdout.flush()

    gourmet.insert_space_charge_markers(kicks_per_line)
    if DEBUG:
        print "after insert_space_charge_markers"
        sys.stdout.flush()
    
    gourmet.complete_setup()
    if DEBUG:
        print "after complete setup"
        sys.stdout.flush()

    if MPI.COMM_WORLD.Get_rank() == 0:
        gourmet.check()
        sys.stdout.flush()


    # check closure for safety
    context = physics_toolkit.BeamlineContext(gourmet.get_initial_particle(),
                                              gourmet.beamline)

    if context.isRing():
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "I am a RING!!!!"
            print "CHEF: Horiz Frac Tune: %.10g" % context.getHorizontalFracTune()
            print "CHEF: Vert Frac Tune: %.10g" % context.getVerticalFracTune()
            print "CHEF: Horiz Frac Eigen Tune: %.10g" % context.getHorizontalEigenTune()
            print "CHEF: Vert Frac Tune: %.10g" % context.getVerticalEigenTune()
            sys.stdout.flush()
    else:
        if MPI.COMM_WORLD.Get_rank() ==0:
            print "NOT a ring!!!  Garbage may result!!!"
            sys.stdout.flush()
    

    (alpha_x, alpha_y, beta_x, beta_y) = synergia.matching.get_alpha_beta(gourmet)
    (tune_x, tune_y, tune_z)           = synergia.matching.get_tunes(gourmet)
    
    if MPI.COMM_WORLD.Get_rank() == 0:
        print "(alpha_x, alpha_y, beta_x, beta_y) = %g, %g, %g, %g" % (alpha_x, alpha_y, beta_x, beta_y)
        print "(tune_x, tune_y, tune_z) = %.10g, %.10g,  %.10g" % (tune_x, tune_y, 1.-tune_z)

    # defining beam_parameters
    charge = 1.0
    initial_phase = 0.0

    beam_parameters = synergia.Beam_parameters(mass,
                                               charge,
                                               kinetic_energy,
                                               initial_phase,
                                               scaling_frequency_Hz=scaling_frequency,
                                               transverse=myopts.get("transverse"),
                                               adjust_zlength_to_freq=myopts.get("adjust_zlength_to_freq"))

    betagamma=beam_parameters.get_beta()*beam_parameters.get_gamma() 

    emittance = myopts.get("emittance")
    pz = beam_parameters.get_gamma() * beam_parameters.get_beta() * beam_parameters.mass_GeV

    (xwidth,xpwidth,rx) = synergia.matching.match_twiss_emittance(emittance,alpha_x,beta_x)
    beam_parameters.x_params(sigma = xwidth, lam = xpwidth * pz,r = rx,offset=xoffset)
    
    (ywidth,ypwidth,ry) = synergia.matching.match_twiss_emittance(emittance,alpha_y,beta_y)
    beam_parameters.y_params(sigma = ywidth, lam = ypwidth * pz,r = ry,offset=yoffset)
    
    zoffset = myopts.get("zoffset")
    zrms=myopts.get("bunchlen")
    zpwidth=myopts.get("dpop")

    # bunch_sp is not used in the current regime
    #bunch_sp=2.0*math.pi*beta*synergia.physics_constants.PH_MKS_c/beam_parameters.get_omega()

    beam_parameters.z_params(sigma = zrms,
                             lam = myopts.get("dpop")* pz,
                             offset = zoffset,
                             offset_p = myopts.get("dpopoffset")*pz)

    
    
    

    #Note! the input term  dpop is in fact (-duop)

    s = 0.0
    line_length = gourmet.orbit_length()
    bunch_spacing = line_length/588.0

    if MPI.COMM_WORLD.Get_rank() ==0:
        print "space_charge =",space_charge
        print "impedance =",impedance

        print "bunch_spacing =",bunch_spacing


    tau = 0.5*line_length/kicks_per_line
    kick_time = 0.0
    beta = beam_parameters.get_beta()
    bunchnp = myopts.get("bunchnp")
    current = bunchnp * \
        synergia.physics_constants.PH_MKS_e/ \
        (bunch_spacing/(beta*synergia.physics_constants.PH_MKS_c))
    if MPI.COMM_WORLD.Get_rank() ==0:
        print "current =",current, " (but not used)"
    
    tgridnum = myopts.get("tgridnum")
    lgridnum = myopts.get("lgridnum")
    griddim = (tgridnum,tgridnum,lgridnum)
    num_particles = int(griddim[0]*griddim[1]*griddim[2] * part_per_cell)

    # creating the bunch
    numbunches = myopts.get("bunches") # I can only deal with one bunch now
    bunches = []
    for bunchnum in range(0,numbunches):
        diag=synergia.Diagnostics(gourmet.get_initial_u(),short=True,save_period=100)
        # if resuming, read the bunch in from the resumedir
        if do_resume:
            bunches.append(s2_fish.Macro_bunch.read_bunch(
                os.path.join(resume_dir, checkpoint_name),
                bunchnum, diagnostics=diag))
        else:
            # otherwise create the new bunch
            bunches.append(s2_fish.Macro_bunch.gaussian(bunchnp, num_particles,
                                                        beam_parameters,
                                                        bucket_num=2*bunchnum,
                                                        diagnostics=diag,
                                                        periodic=myopts.get("periodic")))
        bunches[bunchnum].write_particles("begin-%02d"%bunchnum)

   
    if MPI.COMM_WORLD.Get_rank() ==0:
        print " **********************************************************************" 
        print "beam information:"
        print
        print "emittance = ", emittance
        print "xwidth = ", xwidth
        print "ywidth = ", ywidth
        print "transverse=",beam_parameters.transverse
        if beam_parameters.adjust_zlength_to_freq:
            print "adjusted length =",beam_parameters.z_length
        else:
            print "bunch length (not adjusted)= ", beam_parameters.z_length
           
        print "grid size: ", griddim
        print "number of macroparticles (num_particles): ", num_particles
        print "Charge of bunch: ", bunchnp
        print "Number of bunches: ", numbunches
        
        sys.stdout.flush()


    # solver and geometry

    pipe_radius = myopts.get("pipe_radius")
    space_charge=myopts.get("space_charge")
    if space_charge:
        solver=myopts.get("solver")
        sp_ch_obj = s2_fish.SpaceCharge(solver,griddim,radius_cylindrical=pipe_radius,periodic=myopts.get("periodic"))
        if MPI.COMM_WORLD.Get_rank() ==0:
            print " sp_ch grid=",sp_ch_obj.get_grid()
            print " sp_ch solver=",sp_ch_obj.get_solver()
            print " sp_ch pipe radius=",sp_ch_obj.get_radius_cylindrical()
    else:
        sp_ch_obj=None    	
        
    if MPI.COMM_WORLD.Get_rank() ==0:
        log = open("log", "w")
        output = "start propagation"
        print output
        log.write("%s\n" % output)
        log.flush()

    # if resuming a job, get the starting turn number
    if do_resume:
        rsf = open(os.path.join(resume_dir, "last_turn"))
        turn = int(rsf.readline())
        rsf.close()
    else:
        turn = 0

    # Do diagnostics before propagation
    for bnch in bunches:
        bnch.add_diagnostics(s)
    for jobturn in range(myopts.get("turns")):
        t1 = time.time()
        s = synergia.propagate(s,gourmet,
            bunches, space_charge=sp_ch_obj,
            impedance=None, add_diagnostics=False)
        # Add diagnostics at the end of each turn
        for bnch in bunches:
            bnch.add_diagnostics(s)
        if MPI.COMM_WORLD.Get_rank() ==0:
            output = "turn %d time = %g"%(turn,time.time() - t1)
            print output
            log.write("%s\n" % output)
            log.flush()

        # just completed a turn
        turn = turn+1

        # dump particles if requested (this is different than checkpointing)
        if do_dump and ((turn % dump_freq)==0):
            for bunchnum in range(0,numbunches):
                bunches[bunchnum].write_particles("b-%02d-turn-%05d.h5" %
                                                  (bunchnum,turn))


        # write checkpoint information
        if do_checkpoint and ((turn % checkpoint_freq)==0):
            for bnch in bunches:
                bnch.write_bunch(checkpoint_name)
            # save the last turn checkpointed for resume
            if MPI.COMM_WORLD.Get_rank() == 0:
                ckf = open("last_turn", "w")
                ckf.write("%d\n" % turn)
                ckf.close()


    # save diagnostics after propagation of all turns
    for bunchnum in range(0,numbunches):
        if MPI.COMM_WORLD.Get_rank() == 0:
            bunches[bunchnum].diagnostics.write_hdf5("mi-%02d"%bunchnum)
    for bunchnum in range(0,numbunches):
        bunches[bunchnum].write_particles("end-%02d"%bunchnum)

    if MPI.COMM_WORLD.Get_rank() ==0:
        log.close()

    if MPI.COMM_WORLD.Get_rank() ==0:
        print "elapsed time =",time.time() - t0
