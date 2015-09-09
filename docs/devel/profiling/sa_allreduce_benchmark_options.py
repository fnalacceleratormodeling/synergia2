#!/usr/bin/env python

from synergia_workflow import Options, Job_manager

opts = Options("allreduce_benchmark")
opts.add("gridnum", 32, "size of each of three array dimensions")
opts.add("iterations", 20, "number of iterations")

job_mgr = Job_manager("sa_allreduce_benchmark", opts, standalone=True)

