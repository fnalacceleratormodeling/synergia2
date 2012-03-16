#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("benchmark")
#opts.add("spacecharge", 1, "which space charge calculation: 0=None, 1=Hockney")
opts.add("gridx", 32, "space charge grid cells in x-direction")
opts.add("gridy", 32, "space charge grid cells in y-direction")
opts.add("gridz", 64, "space charge grid cells in z-direction")
opts.add("partpercell", 10, "number of macro particles per grid cell")
opts.add("sortperiod", 1000, "sort period")
opts.add("autotune", True, "automatically tune communication routines")
opts.add("chargecomm", 0, "charge density comm parameter if autotune false, 0 for default")
opts.add("efieldcomm", 0, "electric field comm parameter if autotune false, 0 for default")
opts.add("avoid", True, "whether to use communication avoidance")
opts.add("diagnostics", False, "write diagnostics")
opts.add("verbosity", 1, "verbosity (0=minimal output ... 6=detailed output)")
job_mgr = Job_manager("benchmark", opts,
                      ["cxx_covariance_matrix.xml", "cxx_lattice.xml",
                       "cxx_means.xml"],
                       standalone=True)

