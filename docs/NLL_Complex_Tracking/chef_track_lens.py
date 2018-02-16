#!/usr/bin/env python

import sys, os
import numpy as np
import synergia

# generate test particles for testing nll element.
# 2.5 MeV protons
# normalized emittance 0.3e-6 m-rad
# betax = betay = 0.65

mp = synergia.foundation.pconstants.mp
etot = 0.0025 + mp
normemit = 0.3e-6
betax = 0.65
betay = 0.65

particles = np.loadtxt("particles.txt")
npart = particles.shape[0]

commxx = synergia.utils.Commxx()

refpart = synergia.foundation.Reference_particle(1, mp, etot)
bunch = synergia.bunch.Bunch(refpart, npart, 1.0e10, commxx)
local_particles = bunch.get_local_particles()
local_particles[0:npart,0:4] = particles

knll = 5.479576037e-06
cnll = 0.008105461952

lattice = synergia.lattice.Lattice("channel", synergia.lattice.MadX_adaptor_map())
nllens = synergia.lattice.Lattice_element("nllens", "nl")
nllens.set_double_attribute("knll", knll)
nllens.set_double_attribute("cnll", cnll)
nllens.set_double_attribute("l", 0.0)
lattice.append(nllens)
lattice.set_reference_particle(refpart)

diag = synergia.bunch.Diagnostics_particles("newparticles.h5")

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
#bunch_simulator.add_per_turn(diag)

stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)
propagator = synergia.simulation.Propagator(stepper)

propagator.propagate(bunch_simulator, 1, 1, 1)

np.savetxt("newparticles.txt", local_particles[:,0:4])
