#!/usr/bin/env python

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("foborodobo32_accel")

opts.add("seed", 12345791, "Pseudorandom number generator seed", int)

opts.add("momentum", 2.0, "proton beam momentum")
opts.add("rf_volt", 0.05, "RF voltage [MV]")
opts.add("lag", 0.0, "rf cavity phase in units of 2*pi")

opts.add("matching", "6dmoments", "matching procedure 6dmoments or normal_form")

opts.add("stdx", 1.0e-3, "x standard deviation")
opts.add("stdy", 1.0e-3, "y standard deviation")
opts.add("stdz", 0.055707, "c*dt standard deviation")

opts.add("macro_particles", 100, "number of macro particles")
opts.add("real_particles", 5.0e10, "number of real particles (bunch charge)")

opts.add("turns", 10, "number of turns")

opts.add("spacecharge", None, "space charge [off|2d-openhockney|2d-bassetti-erskine|3d-openhockney", str)
opts.add("gridx", 32, "x grid size")
opts.add("gridy", 32, "y grid size")
opts.add("gridz", 128, "z grid size")

opts.add("xtune", None, "adjust x tune", float)
opts.add("ytune", None, "adjust y tune", float)

opts.add("stepper", "independent", "which stepper to use independent|elements|splitoperator")
opts.add("steps", 128, "# steps")

opts.add("tracks", 100, "number of particles to track")
opts.add("particles", 0, "if non-zero, num particles to save")
opts.add("particles_period", 1, "save  particles every n turns")

job_mgr = synergia_workflow.Job_manager("foborodobo32.py", opts, ["foborodobo32.madx","memusage.py"])
