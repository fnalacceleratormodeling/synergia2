#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("cold_fodo")
opts.add("macro_particles", 10000, "Number of macro particles", int)
opts.add("aperture", "elliptical", "Aperture type (either 'circular', 'elliptical', or 'rectangular'")
opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("real_particles", 1.2e12, "Number of real particles", float)
opts.add("verbose", True, "Verbose propagation", bool)
opts.add("steps", 2, "Number of steps per \"turn\" (really cell)", int)
opts.add("turns", 4, "Number of turns", int)
opts.add("map_order", 2, "Map order", int)
opts.add("emit", 1e-6, "Horizontal and vertical emittance [m rad]", float)
opts.add("stdz", 0.01, "RMS longitudinal length [m]", float)
opts.add("dpop", 1e-4, "delta p / p", float)

job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["fodo.lat"])

