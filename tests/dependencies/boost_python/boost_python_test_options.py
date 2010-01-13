#!/usr/bin/env python

import workflow

opts = workflow.Options("boost_python_test")
opts.add("x", 24.2, "a dummy variable", float)

job_mgr = workflow.Job_manager("boost_python_test.py", opts,["bp_hello.so"])
