#!/usr/bin/env python

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("foborodobo128")

opts.add("seed", 12345791, "Pseudorandom number generator seed", int)

opts.add("matching", "normal_form", "matching procedure 6dmoments or normal_form")

opts.add("stdx", 1.0e-3, "x standard deviation")
opts.add("stdy", 1.0e-3, "y standard deviation")
opts.add("stdz", 10.0, "c*dt standard deviation")

opts.add("macro_particles", 10000, "number of macro particles")
opts.add("real_particles", 5.0e10, "number of real particles (bunch charge)")

opts.add("turns", 10, "number of turns")

opts.add("spacecharge", True, "space charge on or off")
opts.add("gridx", 32, "x grid size")
opts.add("gridy", 32, "y grid size")
opts.add("gridz", 128, "z grid size")

opts.add("stepper", "splitoperator", "which stepper to use independent|elements|splitoperator")
opts.add("steps", 512, "# steps")

opts.add("particles", False, "whether to save particles")
opts.add("tracks", 100, "number of particles to track")

job_mgr = synergia_workflow.Job_manager("foborodobo128.py", opts, ["foborodobo128.madx"])
