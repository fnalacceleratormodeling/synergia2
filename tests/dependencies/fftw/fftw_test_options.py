#!/usr/bin/env python

import workflow

opts = workflow.Options("fftw_test")

job_mgr = workflow.Job_manager("fftw_test.py", opts,["trivial_fftw.so"])
