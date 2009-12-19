#!/usr/bin/env python

import workflow

opts = workflow.Options("fodo")
opts.add("gridnum", 24, "number of grid points to be used for all directions", int)
opts.add("solver", "3d", "solver", str)
opts.add("xoffset", 0.0, "x offset", float)
opts.add("impedance", 0, "whether to use resistive wall kicks", int)
opts.add("piperadius", 0.01, "pipe radius for impedance", float)
opts.add("pipeconduct", 1.4e6,
    "conductivity for pipe [/s], default is for stainless steel", float)
opts.add("spacecharge", 1, "whether to use space charge kicks", int)        
opts.add("np", 2.0e11, "number of particles in real bunch", float)
opts.add("partpercell", 4, "particles per grid cell", float)
opts.add("xwidth", 0.004, "initial horizontal beam width in meters", float)
opts.add("kickspercell", 10, "space-charge kicks per cell", int)
opts.add("dpop", 1.0e-4, "delta p/p", float)

job_mgr = workflow.Job_manager("fodo.py", opts, ["fodo.lat",
                                                 "diagnostics_file.py",
                                                 "fodo-np2e11.h5"])

