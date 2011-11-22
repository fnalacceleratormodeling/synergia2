#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("booster")
opts.add("verbose", False, "Verbose propagation", bool)
opts.add("num_steps", 24*6, "Number of steps per turn", int)
opts.add("num_steps_else", 24*2, "Number of steps per turn", int)
opts.add("num_turns", 2, "Number of turns", int)



opts.add("map_order", 1, "Map order", int)
#opts.add("freq",37866878.4764,"frequency (Hz)",float)
opts.add("num_buckets",84,"number of buckets",int)
#opts.add("harmno",84,"Harmonic number", int)
opts.add("rf_voltage", 0.6/18.0, "RF cavity voltage in MV", float)
opts.add("rfphase",0.0,"rf cavity phase (rad)",float)
#opts.add("latticefile","booster_classic.lat","lattice file",str)
#opts.add("latticefile","booster_0.5.lat","lattice file",str)
opts.add("latticefile","booster_petrenko.lat","lattice file",str)
#opts.add("latticefile","booster_petrenko_leo.lat","lattice file",str)
#opts.add("latticefile","booster_petrenko_zero.lat","lattice file",str)
opts.add("wakefile1","LamD_wake.dat","wake function file",str)
opts.add("wakefile2","LamF_wake.dat","wake function file",str)

opts.add("xrms",0.00856274," RMS x length [m]",float)
opts.add("yrms",0.00329933777917 ," RMS y length [m]",float)  
opts.add("zrms",0.88,"RMS longitudinal length [m]",float)  
opts.add("emit",2.20e-6,"(x-,y-)emittance",float)
#opts.add("norm_emit",5.89533703303356e-07, "Horizontal and vertical emittance [m rad]", float)

opts.add("num_bunches", 1, "number of bunches", int)
opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("num_macro_particles", 100000, "Number of macro particles", int)
opts.add("num_real_particles", 5.e10, "Number of real particles", float)
opts.add("x_offset", 0.00001, "Bunch offset in x", float)
opts.add("y_offset", 0.00001, "Bunch offset in y", float)
opts.add("z_offset", 0.0, "Bunch offset in z", float)


opts.add("impedance", 0, "", int)
opts.add("space_charge", 0, "", int)
opts.add("scgrid",[32,32,64],"Space charge grid",int)
opts.add("grid_ref_distance",0.01,"final spc_grid grid will be scgrid*radius/grid_ref_distance , (m)",float)
opts.add("apertureF",0.021,"average F magnet self-aperture (m)",float)
opts.add("apertureD",0.029,"average D magnet self-aperture (m)",float)
opts.add("apertureL",0.041, "average L pipe self-aperture (m)",float)

opts.add("set_corrections", 1, "", int)
opts.add("correction_file","machine_parameters.sdds",str)

job_mgr = synergia_workflow.Job_manager("booster.py", opts, [opts.get("latticefile"), opts.get("wakefile1"),opts.get("wakefile2"),opts.get("correction_file")], extra_opt_dirs=None)
    