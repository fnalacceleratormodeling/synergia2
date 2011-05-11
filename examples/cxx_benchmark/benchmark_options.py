#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("benchmark")
opts.add("spacecharge", True, "perform space charge calculation")
opts.add("gridx", 32, "space charge grid cells in x-direction")
opts.add("gridy", 32, "space charge grid cells in y-direction")
opts.add("gridz", 64, "space charge grid cells in z-direction")
opts.add("partpercell", 10, "number of macro particles per grid cell")

job_mgr = Job_manager("benchmark", opts,
                      ["cxx_covariance_matrix.xml", "cxx_lattice.xml",
                       "cxx_means.xml"],
                       standalone=True)

