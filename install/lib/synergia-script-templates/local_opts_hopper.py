#!/usr/bin/env python

# The local_options file must be named local_opts.py and placed
# in the Synergia2 job manager search path.

from synergia_workflow import options

# Any instance of the Options class will be added as a suboption.
opts = options.Options('local')
#opts.add("ompnumthreads",2,"number of openmp threads")


# Any instance of the Override class will be used to override
#   the defaults of the existing options.
override = options.Override()
override.numproc = 24
override.procspernode=24
# fix setupsh to point to your own installation
override.setupsh="/global/homes/e/egstern/hopper/synergia2-refactor/setup.sh"
override.template="job_example_hopper"
override.resumetemplate="resume_example_hopper"
override.queue="debug"
override.walltime="00:30:00"
