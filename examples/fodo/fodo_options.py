#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("fodo")
opts.add("num_macro_particles", 320, "Number of macro particles", int)
opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("num_real_particles", 1.2e12, "Number of real particles", float)
opts.add("verbose", True, "Verbose propagation", bool)
opts.add("num_steps", 8, "Number of steps per turn", int)
opts.add("num_turns", 4, "Number of turns", int)
opts.add("map_order", 2, "Map order", int)
opts.add("emit", 1e-6, "Horizontal and vertical emittance [m rad]", float)
opts.add("stdz", 0.01, "RMS longitudinal length [m]", float)
opts.add("dpop", 1e-4, "delta p / p", float)

#opts.add(")
#opts.add("gridnum", 24, "number of grid points to be used for all directions", int)
#opts.add("solver", "3d", "solver", str)
#opts.add("xoffset", 0.0, "x offset", float)
#opts.add("impedance", 0, "whether to use resistive wall kicks", int)
#opts.add("piperadius", 0.01, "pipe radius for impedance", float)
#opts.add("pipeconduct", 1.4e6,
#    "conductivity for pipe [/s], default is for stainless steel", float)
#opts.add("spacecharge", 1, "whether to use space charge kicks", int)
#opts.add("np", 2.0e11, "number of particles in real bunch", float)
#opts.add("partpercell", 4, "particles per grid cell", float)
#opts.add("xwidth", 0.004, "initial horizontal beam width in meters", float)
#opts.add("kickspercell", 10, "space-charge kicks per cell", int)
#opts.add("dpop", 1.0e-4, "delta p/p", float)

job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["fodo.lat",
                                                 "diagnostics_file.py"])

