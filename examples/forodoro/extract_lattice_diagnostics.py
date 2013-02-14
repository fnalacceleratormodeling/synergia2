#!/usr/bin/env python

import synergia

resume = synergia.simulation.Resume()
content = resume.get_content()
lattice_diagnostics = synergia.lattice.Lattice_diagnostics(content.lattice,
                                                            "deposited_charge.h5",
                                                            "deposited_charge")
lattice_diagnostics.update_and_write()
