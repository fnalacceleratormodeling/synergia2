#!/usr/bin/env python
# -*- coding: utf-8 -*-

# The local_options file must be named local_opts.py and placed
# in the Synergia2 job manager search path.

from synergia_workflow import options

# Any instance of the Options class will be added as a suboption.
opts = options.Options('local')
#opts.add("ompnumthreads",2,"number of openmp threads")

# Any instance of the Override class will be used to override
#   the defaults of the existing options.
override = options.Override()
override.template="job_tev"
override.resumetemplate="resumejob_tev"
#override.numprocs = 4

