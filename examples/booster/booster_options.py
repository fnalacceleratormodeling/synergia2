#!/usr/bin/env python
# -*- coding: utf-8 -*-

import synergia_workflow

opts = synergia_workflow.Options("booster")
opts.add("template","job_tev","job template",str)
opts.add("verbosity", 1, "Verbose propagation", int)
opts.add("num_steps", 24*6, "Number of steps per turn", int)
opts.add("num_steps_else", 24*4, "Number of steps per turn", int)
opts.add("num_turns", 2, "Number of turns", int)
opts.add("steps_per_dmag", 1, "Number of steps per D magnet", int)
opts.add("steps_per_fmag", 1, "Number of steps per F magnet", int)


opts.add("map_order", 1, "Map order", int)
#opts.add("freq",37866878.4764,"frequency (Hz)",float)
opts.add("num_buckets",84,"number of buckets",int)
opts.add("full_machine", False, "consider all buckets occupied", bool)
opts.add("wave_number",[0,0,0],"wave number for multibunch instability, trains modulation",int)
opts.add("initial_wave_in_train",False, "if to start with an initial wave in train ", bool)
opts.add("wave_in_train",[0,0,0],"wave number in the train, initial condition",int)
opts.add("zgrid",1000,"longitudinal grid for wake",int)
opts.add("registred_turns",15, "number of previous turns considered for wake",int) 
#opts.add("harmno",84,"Harmonic number", int)
opts.add("rf_voltage", 0.6/18.0, "RF cavity voltage in MV", float)
opts.add("rfphase",0.0,"rf cavity phase (rad)",float)


#opts.add("latticefile","booster_classic.lat","lattice file",str)
#opts.add("latticefile","booster_0.5.lat","lattice file",str)
#opts.add("latticefile","booster_petrenko.lat","lattice file",str)
#opts.add("latticefile","booster_Leo_sbend.lat","lattice file",str)
#opts.add("latticefile","booster_petrenko_leo.lat","lattice file",str)
#opts.add("latticefile","booster_petrenko_zero.lat","lattice file",str)
#opts.add("latticefile","booster_zerochrom.lat","lattice file",str)
#opts.add("latticefile","booster_newtunes.lat","lattice file",str)
opts.add("latticefile","booster_2012.lat","lattice file",str)

#opts.add("wakefileD","BoosterD_g142_muyury.dat","wake function file",str)
#opts.add("wakefileF","BoosterF_g142_muyury.dat","wake function file",str)
opts.add("wakefileD","Dwake.dat","wake function file",str)
opts.add("wakefileF","Fwake.dat","wake function file",str)

#opts.add("wakefileD","BoosterD_g142_mu100.dat","wake function file",str)
#opts.add("wakefileF","BoosterF_g142_mu100.dat","wake function file",str)
#opts.add("wakefileD","LamD_wake.dat","wake function file",str)
#opts.add("wakefileF","LamF_wake.dat","wake function file",str)


opts.add("xrms", 0.005," RMS x length [m]",float) #good ones!!
opts.add("yrms",0.006 ," RMS y length [m]",float)  #good ones!!
opts.add("zrms",0.4 ,"RMS longitudinal length [m]",float) 



#opts.add("emit",2.20e-6,"(x-,y-)emittance",float)
#opts.add("norm_emit",5.89533703303356e-07, "Horizontal and vertical emittance [m rad]", float)

opts.add("num_bunches", 1, "number of bunches", int)
opts.add("seed", 13, "Pseudorandom number generator seed", int)
opts.add("num_macroparticles", 100000, "Number of macro particles", int)
opts.add("num_real_particles", 5e10, "Number of real particles", float)
opts.add("periodic", 1, "periodic bunches along z", int)
opts.add("x_offset", 0.0005, "Bunch offset in x", float)
opts.add("y_offset", 0.0005, "Bunch offset in y", float)
opts.add("z_offset", 0.0, "Bunch offset in z", float)


opts.add("adjust_tunes", 0, "", int)
opts.add("tuneH", 0.78," desired horizontal tune, fractional (0<tuneH<1)",float) 
opts.add("tuneV", 0.85," desired vertical tune, fractional (0<tuneH<1)",float) 

opts.add("adjust_chromaticity", 0, "", int)
opts.add("chromH", -10.," desired horizontal chromaticty",float) 
opts.add("chromV", -5.," desired vertical chromaticty",float) 


opts.add("bpms", 1, "", int)
opts.add("impedance", 0, "", int)
opts.add("space_charge", 0, "", int)
opts.add("aperture",0, "", int)
opts.add("chef_propagate",0, " chef propagate", int)
opts.add("chef_map",0, " chef map", int)
opts.add("turn_period",100, "write particles every turn_period turn", int)
opts.add("scgrid",[32,32,64],"Space charge grid",int)
opts.add("scgrid_L",[128,128,64],"Space charge grid for long section",int)
opts.add("use_comm_divider", 1, " calculate space_charge on all subcomms", int)
opts.add("spc_comm_size", 12, " optimal number of proccs for space_charge, machine dependent", int)
#opts.add("scgrid",[16,16,16],"Space charge grid",int)
#opts.add("scgrid_L",[16,16,16],"Space charge grid for long section",int)
#opts.add("scgrid_L",[512,512,64],"Space charge grid for long section",int)
opts.add("grid_ref_distance",0.01,"final spc_grid grid will be scgrid*radius/grid_ref_distance , (m)",float)
opts.add("apertureF",0.021,"average F magnet self-aperture (m)",float)
opts.add("apertureD",0.029,"average D magnet self-aperture (m)",float)
#opts.add("apertureL",0.041, "average L pipe self-aperture (m)",float)
opts.add("apertureL",0.05715, "average L pipe self-aperture (m)",float)

opts.add("checkpointperiod", 50, "Number of turns to run between checkpoints", int)
opts.add("maxturns",1000, "Maximum number of turns to run before checkpointing and quitting", int)
opts.add("concurrentio", 8, "Maximum number of current io threads for checkpointing", int)
#opts.add("set_corrections", 0, "", int)
#opts.add("correction_file","sextupole_correctors.txt",str)

#job_mgr = synergia_workflow.Job_manager("booster.py", opts, [opts.get("latticefile"), opts.get("wakefileF"),opts.get("wakefileD"),opts.get("correction_file")], extra_opt_dirs=None)

job_mgr = synergia_workflow.Job_manager("booster.py", opts, [opts.get("latticefile"), opts.get("wakefileF"),opts.get("wakefileD")], extra_opt_dirs=None)
      
