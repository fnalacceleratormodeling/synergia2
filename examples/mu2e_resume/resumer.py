#!/usr/bin/env python

import synergia
from resumer_options import opts

resume_dir = "%s/checkpoint" % opts.resume_dir
r = synergia.simulation.Resume(resume_dir)

if opts.max_turns == -1:
    new_max_turns = False
else:
    new_max_turns = True
if opts.verbosity == -1:
    new_verbosity = False
else:
    new_verbosity = True

# these arguments will tell Resume to use the original options
#r.propagate(False, -1, False, -1)
r.propagate(new_max_turns, opts.max_turns, new_verbosity, opts.verbosity)

#lattice_diagnostics = synergia.lattice.Lattice_diagnostics(synergia_lattice,
#                "lattice_deposited_charge.h5", "deposited_charge")
#lattice_diagnostics.update_and_write()

