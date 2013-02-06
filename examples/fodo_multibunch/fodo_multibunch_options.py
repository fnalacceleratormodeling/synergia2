#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("fodo")
opts.add("num_bunches", 3, "number of bunches")
opts.add("bunch_spacing", 1.5, "bunch spacing [m]")
opts.add("macro_particles", 320, "Number of macro particles")
opts.add("seed", 0, "Pseudorandom number generator seed")
opts.add("real_particles", 1.2e12, "Number of real particles")
opts.add("stepper", "independent",
         "Simulation stepper, either 'independent' or 'splitoperator'")
opts.add("step_tracks", 0, "Number of particles to track each step")
opts.add("step_full2", True, "Whether to do full2 diagnostics each step")
opts.add("step_particles", False, "Whether to save all particles each step")
opts.add("turn_tracks", 0, "Number of particles to track each turn")
opts.add("turn_full2", True, "Whether to do full2 diagnostics each turn")
opts.add("turn_particles", False, "Whether to save all particles each turn")
opts.add("verbosity", 2, "Verbosity of propagation")
opts.add("steps", 8, "Number of steps per turn")
opts.add("turns", 4, "Number of turns")
opts.add("checkpointperiod", 2, "Number of turns to run between checkpoints")
opts.add("concurrentio", 8, "Maximum number of current io threads for checkpointing")
opts.add("maxturns", 0, "Maximum number of turns to run before checkpointing and quitting")
opts.add("map_order", 2, "Map order")
opts.add("emit", 1e-6, "Horizontal and vertical emittance [m rad]")
opts.add("stdz", 0.01, "RMS longitudinal length [m]")
opts.add("dpop", 1e-4, "delta p / p")

job_mgr = synergia_workflow.Job_manager("fodo_multibunch.py", opts, ["fodo.lat"])

