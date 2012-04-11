#!/usr/bin/env python

import synergia_workflow

opts = synergia_workflow.Options("mu2e")
opts.add("verbosity", 0, "Verbosity of propagation", int)
opts.add("max_turns", 0, "Maximum number of turns to run before checkpointing and quitting", int)
opts.add("resume_dir", "run.00", "resume directory", str)

job_mgr = synergia_workflow.Job_manager("resumer.py", opts, ["ramp_modules.py"])

