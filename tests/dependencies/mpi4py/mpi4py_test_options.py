#!/usr/bin/env python

import workflow

opts = workflow.Options("mpi4py_test")

job_mgr = workflow.Job_manager("mpi4py_test.py", opts)
