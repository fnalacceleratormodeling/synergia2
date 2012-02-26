#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("fodo")
opts.add("radius", 0.1, "aperture radius [m]", float)
opts.add("macro_particles", 320, "Number of macro particles", int)
opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("real_particles", 1.2e12, "Number of real particles", float)
opts.add("stepper", "independent",
         "Simulation stepper, either 'independent' or 'splitoperator'", str)
opts.add("step_tracks", 0, "Number of particles to track each step", int)
opts.add("step_full2", True, "Whether to do full2 diagnostics each step", bool)
opts.add("step_particles", False, "Whether to save all particles each step", bool)
opts.add("turn_tracks", 0, "Number of particles to track each turn", int)
opts.add("turn_full2", True, "Whether to do full2 diagnostics each turn", bool)
opts.add("turn_particles", False, "Whether to save all particles each turn", bool)
opts.add("verbosity", 0, "Verbosity of propagation", int)
opts.add("steps", 8, "Number of steps per turn", int)
opts.add("turns", 4, "Number of turns", int)
opts.add("max_turns", 0, "Maximum number of turns to run before checkpointing and quitting", int)
opts.add("map_order", 2, "Map order", int)
opts.add("emit", 1e-6, "Horizontal and vertical emittance [m rad]", float)
opts.add("stdz", 0.01, "RMS longitudinal length [m]", float)
opts.add("dpop", 1e-4, "delta p / p", float)

job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["fodo.lat"])

