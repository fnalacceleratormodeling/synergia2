#!/usr/bin/env python

import workflow

opts = workflow.Options("python_test")
opts.add("x", 24.2, "a dummy variable", float)

job_mgr = workflow.Job_manager("python_test.py", opts)
