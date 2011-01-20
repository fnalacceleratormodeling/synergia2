#!/usr/bin/env python

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("mi")
opts.add("num_macro_particles", 1000, "Number of macro particles", int)
opts.add("seed", 0, "Pseudorandom number generator seed", int)
opts.add("num_real_particles", 1.1e11, "Number of real particles", float)
opts.add("verbose", True, "Verbose propagation", bool)
opts.add("num_steps", 728, "Number of steps per turn", int)
opts.add("num_turns", 20, "Number of turns", int)
opts.add("map_order", 1, "Map order", int)
# 95% emittance is 18 mm-mr is geometric emittance*6*pi
opts.add("norm_emit",18.0e-6/(6.0*pi), "Horizontal and vertical emittance [m rad]", float)
opts.add("stdz", 0.3, "RMS longitudinal length [m]", float)
#  RF voltage is 1MV divided among 18 cavities
opts.add("rf_voltage", 1.0/18, "RF cavity voltage in MV", float)
opts.add("x_offset", 0.0, "Bunch offset in x", float)
opts.add("y_offset", 0.0, "Bunch offset in y", float)
opts.add("z_offset", 0.0, "Bunch offset in z", float)

job_mgr = synergia_workflow.Job_manager("mi.py", opts, ["mi20-egs-thinrf.lat"])
