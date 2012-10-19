#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("fodo")
opts.add("macro_particles", 320, "Number of macro particles", int)
opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("real_particles", 1.2e12, "Number of real particles", float)
opts.add("verbose", False, "Verbose propagation", bool)
opts.add("steps", 1, "Number of steps per \"turn\" (really cell)", int)
opts.add("turns", 1, "Number of turns", int)
opts.add("checkpointperiod", 1, "Number of turns to run between checkpoints", int)
opts.add("map_order", 2, "Map order", int)
opts.add("emit", 1e-6, "Horizontal and vertical emittance [m rad]", float)
opts.add("stdz", 0.01, "RMS longitudinal length [m]", float)
opts.add("dpop", 1e-4, "delta p / p", float)

job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["fodo.lat"])

