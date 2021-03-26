#!/usr/bin/env python3

from math import pi
import synergia_workflow

opts = synergia_workflow.Options("channel")

opts.add("order", 7, "normal form order")
opts.add("turns", 1000, "number of turns")
opts.add("skew", 0.0, "skew quad strength (try 0.005 for fun)")
opts.add("sext", 0.0, "strength of sextupole")
opts.add("sksext", 0.0, "strength of skew sextupole")
opts.add("octo", 0.0, "strength of octopole (try 0.05 for fun)")
opts.add("skocto", 0.0, "strength of skew octopole")

opts.add("RFVolt", 1.0, "RF Cavity voltage [MV] (try 0.005 or 10 for fun)")
opts.add("offset", 0.001, "offset for test particles")

job_mgr = synergia_workflow.Job_manager("channel.py", opts, ["channel_template.seq", "ramp_module.py"])
