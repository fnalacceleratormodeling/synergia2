#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("elens")
#opts.add("spacecharge", 1, "which space charge calculation: 0=None, 1=Hockney")
opts.add("lattice_file", "el_Xoforodo.lsx", "lattice file to run")

#opts.add("rf_voltage", 0.00996169973931, "rf voltage") this is for sync tune 1 1/96 with 12 RF cavities

# this is for sync tune 1/13 with 12 rf cavities
opts.add("rf_voltage", 0.261977646285, "rf voltage")

# sync tune 1/13 with 6 RF cavities
#opts.add("rf_voltage", 0.52395529257, "rf voltage")

# this is the product of current*length to counteract the space charge kick
#elens_JL = 23.439596537
elens_JL  = 24.0
opts.add("elensradius", 0.0041525, "electron lens beam radius [m]")
opts.add("elensenergy", 0.010, "electron lens beam energy [MV]")
opts.add("elenslength", 2.0, "electron lens length [m]")
#opts.add("elenscurrent", elens_JL/opts.elenslength/6.0, "electron lens current [A]")
opts.add("elenscurrent", 0.0, "electron lens current [A]")
opts.add("elenslongrms", -1.0, "elens longitudinal RMS pulse length")
opts.add("elensdivide", 1, "divide number of elenses by this factor")

opts.add("harmon", 50, "harmonic number")
opts.add("transport", "maps", "propagation method [maps|chef|libff]")

opts.add("bpms", True, "use bpm")

opts.add("stepper", "splitoperator", "stepper to use: [independent|elements|splitoperator|splitoperatorelements|choice]")

opts.add("space_charge_3dh", 1, "use space charge 3d")
opts.add("space_charge_2dh", 0, "use space charge 2d")
opts.add("space_charge_rec", 0, "use rectangular space charge")

opts.add("solver", "3doh", "space charge solver to use [3doh, 2doh, rect, dummy]")

opts.add("spc_comm_size", 32, "size of space charge communcator")

opts.add("gridx", 64, "space charge grid cells in x-direction")
opts.add("gridy", 64, "space charge grid cells in y-direction")
opts.add("gridz", 128, "space charge grid cells in z-direction")

opts.add("steps_per_quad", 1, "steps per quad")
opts.add("num_steps_else", 100, "number steps other elements")
opts.add("num_steps", 85, "number of steps for independent stepper")
opts.add("map_order", 1, "map order")

#opts.add("xrms", 0.0041525, "x rms")
#opts.add("yrms", 0.00414, "y rms")
opts.add("xrms", 0.0041525, "x rms")
opts.add("yrms", 0.0041525, "y rms")
opts.add("zrms", 0.5, "z rms")

# these lattice functions come from the unperturbed lattice.  Use them to
# calculate the appropriate emittance to use for populating bunches for the
# perturbed lattice so the emittance starts the same.
#beta_x = 17.3259852015
#beta_y = 17.2246414528

# These are the beta functions for equal betax and betay where the
# tune is set to 0.74
beta_x = 17.2336287
beta_y = 17.2336287

opts.add("emitx", opts.xrms**2/beta_x, "x emittance")
opts.add("emity", opts.yrms**2/beta_y, "y emittance")

opts.add("num_turns", 1000, "number of turns")
opts.add("maxturns", 3000, "maximum number of turns in one run")

opts.add("num_bunches", 1, "number of bunches")
opts.add("macroparticles", 100000, "number of macroparticles")
opts.add("real_particles", 0.666e11, "number of charges")

opts.add("load_bunch", 0, "load the bunch")
opts.add("save_bunch", 0, "save the bunch")

opts.add("seed", 13, "random seed")

opts.add("turn_period", 25, "save particles every n turns")
opts.add("turn_track", 10000, "# tracks to save, 0: don't save any")
opts.add("spc_tuneshift", 1, "save tune shift calculation diagnostics")
opts.add("cell_particles", False, "Add diagnostics per cell")

opts.add("xoffset", 0.0, "offset bunch in x")
opts.add("yoffset", 0.0, "offset bunch in y")
opts.add("zoffset", 0.0, "offset bunch in z")

opts.add("errsize", 0.0, "size of introduced error")
opts.add("errelement", 12, "the element that receives the error")

opts.add("checkpointperiod", 50, "save checkpoint every n turns")
opts.add("concurrentio", 8, "use this many processors concurrently for checkpointing")

opts.add("zparticle", False, "particle at 0,0")
opts.add("set_nux", 0.74, "value to set nux for the lattice")
opts.add("set_nuy", 0.74, "value to set nuy for the lattice")

opts.add("verbosity", 1, "verbosity (0=minimal output ... 6=detailed output)")
job_mgr = Job_manager("cxx_elens", opts,
                      ["el_Xoforodo.lsx"],
                       standalone=True)

