#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("fodo")
opts.add("map_order", 1, "map order")
opts.add("steps_per_element", 2, "steps per element")
opts.add("x_emit", 1.0e-6, "x RMS emittance [m-rad]")
opts.add("y_emit", 1.0e-6, "y RMS emittance [m-rad]")
opts.add("z_std", 0.01, "z RMS length [m]")
opts.add("dpop", 1.0e-4, "(delta p)/p")
opts.add("real_particles", 1.2e12, "number of physical particles in bunch")
opts.add("macro_particles", 100, "number of simulation particles")
opts.add("seed", 1415926, "random number seed; 0 for automatic calculation (GSL)")
opts.add("turns", 4, "number of times to track through fodo lattice")
opts.add("max_turns", 4, "maximum number of turns to run before checkpointing and stopping; 0 to not stop")
opts.add("verbosity", 2, "verbosity level of simulation")

# Create the job manager for the simulation fodo.py, including the above options.
# When creating job directories, include the file fodo.lat.
job_mgr = synergia_workflow.Job_manager("fodo.py", opts, ["fodo.lat"])

