#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("offdiag")
#opts.add("spacecharge", 1, "which space charge calculation: 0=None, 1=Hockney")
opts.add("lattice_file", "el_Xoffdiag.lsx", "lattice file to run")

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
opts.add("elenslength", 1.0, "electron lens length [m]")
#opts.add("elenscurrent", elens_JL/opts.elenslength/6.0, "electron lens current [A]")
opts.add("elenscurrent", 0.0, "electron lens current [A]")
opts.add("elenslongrms", 0.5, "elens longitudinal RMS pulse length")
opts.add("elensdivide", 2, "divide number of elenses by this factor")
opts.add("elensadaptive", False, "adaptive elens radius")

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

opts.add("magiccomp", 0, "magic compensation, 0: no compensation, 1: exactly opposite space charge kick, 2: double compensation every other kick")
opts.add("sccomp", 1.0, "space charge compensation compensation factor")

opts.add("steps_per_quad", 1, "steps per quad")
opts.add("num_steps_else", 100, "number steps other elements")
opts.add("num_steps", 72, "number of steps for independent stepper")
opts.add("map_order", 1, "map order")

#opts.add("xrms", 0.0041525, "x rms")
#opts.add("yrms", 0.00414, "y rms")
opts.add("xrms", 0.0041525, "x rms")
opts.add("yrms", 0.0041525, "y rms")
opts.add("zrms", 0.5, "z rms")
opts.add("dpoprms", 1.0e-3, "dp/p rms")

opts.add("k2l", 0.0, "k2l of sextupole element")

# these lattice functions come from the unperturbed lattice.  Use them to
# calculate the appropriate emittance to use for populating bunches for the
# perturbed lattice so the emittance starts the same.
#beta_x = 17.3259852015
#beta_y = 17.2246414528

# These are the beta functions lattice with tunes 0.74 0.84
beta_x = 17.2742174685
beta_y = 17.2779573169

opts.add("beamparams", False, "Use emittance, alpha, beta beam params")
opts.add("transhemit", 1.09607333934e-06, "horizontal emittance for beam")
opts.add("transvemit", 1.04757819876e-06, "vertical emittance for beam")
opts.add("hbeta", 18.3185523495, "horizontal beta for beam")
opts.add("vbeta", 18.0230840973, "vertical beta for beam")
opts.add("halpha", 1.83173893189, "horizontal alpha for beam")
opts.add("valpha", -1.85824588289, "vertical alpha for beam")

opts.add("flatbucket", False, "Generate a bucket that has a flat longitudinal profile")

#opts.add("emitx", opts.xrms**2/beta_x, "x emittance")
#opts.add("emity", opts.yrms**2/beta_y, "y emittance")
opts.add("emitx", 1.00055865e-06, "x emittance")
opts.add("emity", 1.00055865e-06, "y emittance")

opts.add("apertures", False, "Put hard cut on particles")
opts.add("aperture_sigma", 4.0, "Aperture radius in sigmas")

opts.add("num_turns", 1000, "number of turns")
opts.add("maxturns", 3000, "maximum number of turns in one run")

opts.add("num_bunches", 1, "number of bunches")
opts.add("macroparticles", 1000000, "number of macroparticles")
opts.add("real_particles", 2.0e11, "number of charges")

opts.add("load_bunch", 0, "load the bunch")
opts.add("save_bunch", 0, "save the bunch")

opts.add("seed", 13, "random seed")

opts.add("turn_particles", 1, "save particles every turn")
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
opts.add("set_nux", 0.0, "value to set nux for the lattice")
opts.add("set_nuy", 0.0, "value to set nuy for the lattice")

opts.add("verbosity", 1, "verbosity (0=minimal output ... 6=detailed output)")
job_mgr = Job_manager("cxx_offdiag", opts,
                      ["el_Xoffdiag.lsx"],
                       standalone=True)

