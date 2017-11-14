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

#particles = np.loadtxt("particles.txt")
#npart = particles.shape[0]
npart = 2

commxx = synergia.utils.Commxx()

refpart = synergia.foundation.Reference_particle(1, mp, etot)
bunch = synergia.bunch.Bunch(refpart, npart, 1.0e10, commxx)

local_particles = bunch.get_local_particles()
local_particles[1, 0] = .001
local_particles[1, 2] = -0.002

knll = 5.479576037e-06
cnll = 0.008105461952

total_length = 1.0e-6

lattice = synergia.lattice.Lattice("channel", synergia.lattice.MadX_adaptor_map())

d1 = synergia.lattice.Lattice_element("drift", "d1")
d1.set_double_attribute("l", total_length/2.0)
lattice.append(d1)

nllens = synergia.lattice.Lattice_element("nllens", "nl")
nllens.set_double_attribute("knll", knll)
nllens.set_double_attribute("cnll", cnll)
nllens.set_double_attribute("l", 0.0)
lattice.append(nllens)

d2 = synergia.lattice.Lattice_element("drift", "d2")
d2.set_double_attribute("l", total_length/2.0)
lattice.append(d1)

lattice.set_reference_particle(refpart)

#diag = synergia.bunch.Diagnostics_particles("newparticles.h5")

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
#bunch_simulator.add_per_turn(diag)

stepper = synergia.simulation.Independent_stepper(lattice, 1, 1)
propagator = synergia.simulation.Propagator(stepper)

propagator.propagate(bunch_simulator, 1, 1, 1)

print np.array2string(local_particles[0, 0:4])
print np.array2string(local_particles[1, 0:4])

