#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("sis18_frs")

opts.add("comm_divide", 32, "size of communicator")
opts.add("generate_bunch", True, "Generate matched bunch internally", bool)
opts.add("particles_file", "sis18-6_particles.txt", "text file containing bunch particles", str)

# these values came from Frank Schmidt from a mad X run 2012-10-16
opts.add("kqf", 3.10827e-01, "focussing quad k strength", float)
opts.add("kqd",  -4.95996e-01, "defocussing quad k strength", float)

opts.add("k2l", 0.0, "sextupole setting", float)

opts.add("xtune", 0.338, "setting for xtune", float)
opts.add("ytune", 0.200, "setting for ytune", float)
opts.add("tune_tolerance", 1.0e-7, "Tolerance for tune adjustment", float)
opts.add("xml_save_lattice", False, "Whether to save cooked lattice as xml file", bool)
opts.add("emitx", 12.57e-6, "real 2 sigma Horizontal emittance [m rad]", float)
opts.add("emity", 9.3e-6, "real 2 sigma Vertical emittance [m rad]", float)

# dpop of 1.038*2.4e-4/3 gives stdz too large.
#opts.add("dpop", 1.038*2.4e-4/3.0, "Delta-p/p spread", float)
opts.add("dpop", 2.4e-4/3.0, "Delta-p/p spread", float)

opts.add("radius", 0.5, "aperture radius [m]", float)
opts.add("macro_particles", 10485760, "Number of macro particles", int)
opts.add("seed", 4, "Pseudorandom number generator seed", int)
opts.add("real_particles", 2.94e10, "number of real particles", float)
#opts.add("real_particles", 30486095439.0, "Number of real particles", float)
#opts.add("real_particles", 19294997113.0 , "Number of real particles", float)
#opts.add("real_particles", 18733006906.0 , "Number of real particles", float)
opts.add("spacecharge", True, "whether space charge is on", bool)
opts.add("solver", "2dopen-hockney", "solver to use, '2dopen-hockney','3dopen-hockney', '2dbassetti-erskine'", str)

opts.add("matching", "normalform", "Matching procedure: normalform|6dlinear", str)

opts.add("test_particles", True, "include test particles in bunch", bool)
# second particle 0.1 sigma x offset, sigma=.005075
opts.add("x_test", 6.3361e-3/10.0, "x offset of x test particle", float)
# third particle 0.1 sigma y offset sigma = 004497
opts.add("y_test", 5.59708e-3/10.0, "y offset of y test particle", float)
opts.add("test_step", 1.0, "offset in x_test or y_test  between test particles", float)

opts.add("z_offset", 0.0, "z offset of test particles", float)

#========= grid size for SC solver ===========
#opts.add("gridx", 32, "grid points in x for solver", int)
#opts.add("gridy", 32, "grid points in y for solver", int)
#opts.add("gridz", 32, "grid points in z for solver", int)
opts.add("gridx", 128, "grid points in x for solver", int)
opts.add("gridy", 128, "grid points in y for solver", int)
opts.add("gridz", 128, "grid points in z for solver", int)

opts.add("use_maps", "none", "use maps for propagation either all, none, onlyrf, nonrf")
opts.add("allmaps", False, "Use all maps for propagation", bool)
opts.add("stepper", "splitoperator",
         "Simulation stepper, either 'independent','elements','splitoperator','soelements'", str)

#opts.add("scratch", "/scratch", "directory for temporary diagnostic files", str)
opts.add("scratch", "", "directory for temporary diagnostic files", str)
opts.add("turn_tracks", 201, "Number of particles to track each turn", int)
opts.add("turn_full2", True, "Whether to do full2 diagnostics each turn", bool)
opts.add("turn_particles", False, "Whether to save all particles each turn", bool)
opts.add("turn_particles_period", 100, "period for dumping particles", int)

opts.add("step_tracks", 0, "number of particles to track each step", int)
opts.add("step_full2", False, "Whether to do full2 diagnostics each step", bool)
opts.add("step_particles", False, "Whether to save all particles each step", bool)

opts.add("per_operator_full2", False, "full2 diagnostics per operator")
opts.add("per_operator_tracks", 0, "output tracks per operator")
opts.add("per_operator_particles", 0, "output particles per operator")

opts.add("verbosity", 1, "Verbosity of propagation", int)
opts.add("steps", 71, "Number of steps per turn", int)
opts.add("turns", 1024, "Number of turns", int)
opts.add("checkpointperiod", 2000, "Number of turns to run between checkpoints", int)
opts.add("concurrent_io", 32, "number of concurrent io threads for checkpointing", int)
opts.add("maxturns", 15000, "Maximum number of turns to run before checkpointing and quitting", int)
opts.add("map_order", 1, "Map order", int)

opts.add("dump_efield", False, "Whether do dump efield every turn", bool)

opts.add("nsigma", 8.0, "nsigma for solver", float)
opts.add("long_kicks", True, "use longitudinal kicks", bool)
opts.add("cutoffnsigma", 2.6, "cut off bunch at cutoffnsigma standard deviations")
job_mgr = synergia_workflow.Job_manager("sis18_frs.py", opts, ["sis18-6.mad","ramp_module.py"])
