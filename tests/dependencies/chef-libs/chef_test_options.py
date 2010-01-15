#!/usr/bin/env python

import workflow

opts = workflow.Options("chef_test")

job_mgr = workflow.Job_manager("chef_test.py", opts)
