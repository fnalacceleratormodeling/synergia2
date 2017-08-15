#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("space_charge_drift")

opts.add("current", 14.0e-3, "Current [A]")
opts.add("ke", 2.5e-3, "Kinetic Energy [GeV]")
opts.add("nemit", 0.3e-6, "normalized emittance m-rad")
opts.add("betax", 0.5, "effective beta_x")

opts.add("betay", 0.5, "effective beta_y")

opts.add("blen", 0.02, "bunch length [m]")
opts.add("driftlength", 0.1, "drift length [m]")

opts.add("turns", 100, "number of turns")
opts.add("steps", 25, "steps per turn")

opts.add("macroparticles", 409600, "macro particles")

opts.add("seed", 12345679, "seed for random number generation")

opts.add("gridx", 64, "x grid points")
opts.add("gridy", 64, "y grid points")
opts.add("gridz", 32, "z grid points")

opts.add("solver", "2d-kv", "solver to use: 2d-openhockney|2d-kv")

opts.add("centered", False, "is field calculation assumed to be centered on beamline")
opts.add("plot", False, "plot comparison", bool)
opts.add("particles", False, "Whether to save particles")
opts.add("particles_period", 1, "save particles every n turns")
opts.add("verbosity", 1, "chattiness of simulation")

# Create the job manager for the simulation fodo_workflow.py, including the 
# above options. When creating job directories, include the file fodo.lat.
job_mgr = Job_manager("space_charge_drift.py", opts, ["space_charge_drift_options.py"])

