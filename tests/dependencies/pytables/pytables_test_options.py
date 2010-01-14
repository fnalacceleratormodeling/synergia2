#!/usr/bin/env python

import workflow

opts = workflow.Options("pytables_test")

job_mgr = workflow.Job_manager("pytables_test.py", opts)
