#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("fodo_space_charge")
opts.add("map_order", 1, "map order")
opts.add("stepper", "element", 
  "split operator stepper, either fixed steps per element or fixed length steps",
  valid_values=["element", "fixed"])
opts.add("steps_per_element", 2, "steps per element (only for element stepper)")
opts.add("num_steps", 8, "total number of steps (only for fixed length stepper)")
opts.add("space_charge", "2d", "space charge model", 
  valid_values=["none", "bassetti", "2d", "3d"])
opts.add("gridx", 16, "grid cells in x-direction (2d and 3d space charge)")
opts.add("gridy", 16, "grid cells in y-direction (2d and 3d space charge)")
opts.add("gridz", 16, "grid cells in z-direction (2d and 3d space charge)")
opts.add("x_emit", 1.0e-6, "x RMS emittance [m-rad]")
opts.add("y_emit", 1.0e-6, "y RMS emittance [m-rad]")
opts.add("z_std", 0.01, "z RMS length [m]")
opts.add("dpop", 1.0e-4, "(delta p)/p")
opts.add("real_particles", 1.2e12, "number of physical particles in bunch")
opts.add("macro_particles", 50000, "number of simulation particles")
opts.add("seed", 1415926, 
         "random number seed; 0 for automatic calculation")
opts.add("turns", 1, "number of times to track through fodo lattice")
opts.add("max_turns", 0, 
         "maximum number of turns to run before checkpointing and stopping; 0 to not stop")
opts.add("verbosity", 2, "simulation verbosity level")

# Create the job manager for the simulation fodo_space_charge.py, including the 
# above options. When creating job directories, include the file fodo.lat.
job_mgr = Job_manager("fodo_space_charge.py", opts, ["fodo.madx"])

