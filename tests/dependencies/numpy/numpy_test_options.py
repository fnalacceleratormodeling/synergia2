#!/usr/bin/env python

import workflow

opts = workflow.Options("numpy_test")

job_mgr = workflow.Job_manager("numpy_test.py", opts)
