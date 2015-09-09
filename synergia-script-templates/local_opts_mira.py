#!/usr/bin/env python

# The local_options file must be named local_opts.py and placed
# in the Synergia2 job manager search path.

from synergia_workflow import options

# Any instance of the Options class will be added as a suboption.
opts = options.Options('local')

# Any instance of the Override class will be used to override
#   the defaults of the existing options.
override = options.Override()
# change the next line with your ALCF account
override.account = "xxxxx"
override.numproc = 4096
override.procspernode=4
override.template="job_mira"
override.resumetemplate="resume_mira"
override.multitemplate="mira_multijob"
override.resumemultitemplate="mira_resumemulti"
override.queue="default"
override.walltime="06:00:00"
