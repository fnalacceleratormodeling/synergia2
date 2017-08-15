#!/usr/bin/env python
# -*- coding: utf-8 -*-

import synergia_workflow

opts = synergia_workflow.Options("iota")

#lattice, lattice simulator and stepper
opts.add("nnl", 0, "nonlinear insert", int)
opts.add("map_order", 1, "Map order", int)
opts.add("chef_propagate",1, " chef propagate", int)
opts.add("chef_map",0, " chef map", int)
opts.add("steps_per_sbend", 4, "Number of steps per sbend magnets", int)
opts.add("num_steps_straight",288,"Number of steps per straight section", int)
opts.add("num_steps", 1640, "Number of steps for independent stepper", int)
opts.add("tunes_and_chroms",1, "calculates and print tunes and chromaticities", int)
opts.add("harmon",1,"Harmonic number", int)
#opts.add("rf_voltage", 0.6/18.0, "RF cavity voltage in MV", float)
#opts.add("adjust_tunes", 0, "", int)
#opts.add("tuneH", 0.78," desired horizontal tune, fractional (0<tuneH<1)",float) 
#opts.add("tuneV", 0.85," desired vertical tune, fractional (0<tuneH<1)",float) 
#opts.add("adjust_chromaticity", 0, "", int)
#opts.add("chromH", -10.," desired horizontal chromaticty",float) 
#opts.add("chromV", -5.," desired vertical chromaticty",float) 


#apertures
opts.add("if_aperture",1, "bool option (0,1)", int)
opts.add("aperture_bending",0.025, "rectangular half width aperture in bending magnets(m)",float)
opts.add("aperture_straight",0.02375, "radius circular aperture in straight section(m)",float)

#space charge
opts.add("space_charge_rec", 0, "bool option (0,1), space charge rectangular solver", int)
opts.add("space_charge_3dh", 0, "bool option (0,1),space charge 3d Hockney solver", int)
opts.add("space_charge_2dh", 0, "bool option (0,1), space charge 2d Hockney solver", int)
opts.add("scgrid_straight",[64,64,128],"Space charge grid for straight section",int)
opts.add("scgrid_bending",[64,64,128],"Space charge grid for bending section",int)
opts.add("scgrid_3dh",[64,64,128],"Space charge grid for 3d Hockney solver",int)
opts.add("spc_comm_size",32," optimal size of the comunicator used by FFT in the 3d Hockney solver",int)




#impedance
opts.add("impedance", 0, "bool option (0,1), impedance", int)
opts.add("zgrid",1000,"longitudinal grid for wake",int)
opts.add("wakefile_straight","IOTA_straight_rw_wake.dat","wake function file",str)
opts.add("waketype","XLYLZ","wake type",str)
opts.add("full_machine", False, "consider all buckets occupied", bool)
opts.add("wave_number",[0,0,0],"wave number for multibunch instability, trains modulation",int)
opts.add("registred_turns",100, "number of previous turns considered for wake",int)



#bunch
opts.add("xrms", 0.0018," RMS x length [m]",float) 
opts.add("yrms",0.0018 ," RMS y length [m]",float)  
opts.add("zrms",1.7 ,"RMS longitudinal length [m]",float) 
opts.add("x_offset", 0.0006, "Bunch offset in x", float)
opts.add("y_offset", 0.0003, "Bunch offset in y", float)
opts.add("z_offset", 0.0, "Bunch offset in z", float)

opts.add("num_bunches", 1, "number of bunches", int)
opts.add("num_macroparticles", 100000, "Number of macro particles", int)
opts.add("num_real_particles", 1e11, "Number of real particles", float)
opts.add("coasting_beam", 1, "bool option (0,1) for generating a coasting beam",int)
opts.add("load_bunch",0," read bunch for file",int)
opts.add("input_bunch_h5file","iota_input_particles.h5","h5 format file name with input bunch, (format as created by diagnostics_particle) ",str)


opts.add("seed", 13, "Pseudorandom number generator seed", int)




#diagnostics properties

opts.add("turn_period",100, "write particles every turn_period turn", int)
opts.add("bulk_track",0, "bool option (0,1) to track and save particles in bulk_track",int)
opts.add("spc_tuneshift",1, "bool option (0,1) to create a diagnostics for the incoherent tune shift" ,int)
opts.add("apertures_loss",1, "bool option (0,1) to create a loss diagnostics" ,int)

#propagator properties
opts.add("num_turns", 1, "Number of turns", int)
opts.add("checkpointperiod", 50, "Number of turns to run between checkpoints", int)
opts.add("maxturns",3000, "Maximum number of turns to run before checkpointing and quitting", int)
opts.add("concurrentio", 8, "Maximum number of current io threads for checkpointing", int)
opts.add("verbosity", 1, "Verbose propagation", int)


job_mgr = synergia_workflow.Job_manager("iota.py", opts, [opts.get("wakefile_straight")], extra_opt_dirs=None)
      
