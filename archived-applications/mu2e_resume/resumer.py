#!/usr/bin/env python

import sys
import os
import mpi4py.MPI as MPI
import synergia
from resumer_options import opts

myrank = synergia.utils.Commxx().get_rank()

t0 = MPI.Wtime()

resume_dir = "%s/checkpoint" % opts.resume_dir
resume = synergia.simulation.Resume(resume_dir)
content = resume.get_content()

if opts.max_turns == -1:
    new_max_turns = False
else:
    new_max_turns = True
if opts.verbosity == -1:
    new_verbosity = False
else:
    new_verbosity = True

# these arguments will tell Resume to use the original options
#resume.propagate(False, -1, False, -1)
resume.propagate(new_max_turns, opts.max_turns, new_verbosity, opts.verbosity)

t1 = MPI.Wtime()

if myrank == 0:
    print
    print "Propagate time =", t1 - t0

lattice_diagnostics = synergia.lattice.Lattice_diagnostics(content.lattice,
                "lattice_deposited_charge.h5", "deposited_charge")
lattice_diagnostics.update_and_write()

