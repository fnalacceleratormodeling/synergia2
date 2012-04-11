#!/usr/bin/env python

# The local_options file must be named local_opts.py and placed
# in the Synergia2 job manager search path.

from synergia_workflow import options

# Any instance of the Options class will be added as a suboption.
opts = options.Options('local')
#opts.add("ompnumthreads",2,"number of openmp threads")
opts.add("loadbalance",True,"whether to specify --loadbalance for mpirun", bool)


# Any instance of the Override class will be used to override
#   the defaults of the existing options.
override = options.Override()
override.account = "mu2e"
override.numproc = 32
override.procspernode=32
override.setupsh="/home/cspark/synergia2-checkpoint/setup.sh"
override.template="job_example_wilson"
override.queue="intel12"
#override.queue="amd32"
override.walltime="24:00:00"
