#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("fodo_cxx")
opts.add("gridx", 32, "space charge grid cells in x-direction")
opts.add("gridy", 32, "space charge grid cells in y-direction")
opts.add("gridz", 128, "space charge grid cells in z-direction")
opts.add("macroparticles", 1048576, "number of macroparticles")
opts.add("real_particles", 2.94e12, "total charge")
opts.add("turns", 10, "number of turns")

job_mgr = Job_manager("fodo_cxx", opts, [], standalone=True)

