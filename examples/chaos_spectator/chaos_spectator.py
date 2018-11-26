#!/usr/bin/env python

import sys, os
import numpy as np
import matplotlib.pyplot as plt
import synergia

import tables

# The first part of this script reads a lattice containing a FODO cell with
# a sextupole.  It propagates a bunch of particles.  Imagine that you would
# the stability of specific particles.  You create spectator particles
# with coordinates from the specific particle and one displaced by an
# infinitesimal offset.  Run 1000 turns writing out diagnostics.  The second
# part of the script reads the diagnostic files and makes plots.

lattice = synergia.lattice.MadX_reader().get_lattice("model", "henon.madx")

print lattice.as_string()

map_order = 1

real_particles = 1.0e10
macro_particles = 10

commxx = synergia.utils.Commxx()

bunch = synergia.bunch.Bunch(lattice.get_reference_particle(),
                             macro_particles, 6, real_particles, commxx)


lp = bunch.get_local_particles()
print "lp.shape: ", lp.shape
lsp = bunch.get_local_spectator_particles()

for i in range(lp.shape[0]):
    lp[i, 0] = 0.001*i

# spectator particles, lp 1, lp 4, lp 9 with small offsets also

lsp[0, 0:6] = lp[1, 0:6]

lsp[1, 0:6] = lp[1, 0:6]
lsp[1, 1] = 4.0e-16

lsp[2, 0:6] = lp[4, 0:6]
lsp[3, 0:6] = lp[4, 0:6]
lsp[3, 1] = 4.0e-16

lsp[4, 0:6] = lp[9, 0:6]
lsp[5, 0:6] = lp[9, 0:6]
lsp[5, 1] = 4.0e-16

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)

bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_full2("diagnostics.h5"))
bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_particles("particles.h5"), 10)
bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("tracks.h5", 10))
bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_spectator_track("spectracks.h5", 6))

stepper = synergia.simulation.Independent_stepper(lattice, map_order, 1)

propagator = synergia.simulation.Propagator(stepper)

propagator.propagate(bunch_simulator, 1000, 0, 1)

del(propagator)
del(stepper)
del(bunch_simulator)
del(bunch)
del(lattice)

# read and analyze spectator particle diagnostics

def tdiff(c1, c2):
    return np.sqrt(np.dot((c1-c2), (c1-c2)))

def le(t1, t2):
    # t1, t2 are two tracks nx6
    w0 = tdiff(t1[0,:], t2[0, :])
    n = t1.shape[0]
    letrk = np.zeros((n), dtype='d')
    for i in range(n):
        letrk[i] = np.log(tdiff(t1[i,:], t2[i,:])/w0)
    return letrk

h5 = tables.open_file("spectracks.h5")

tracks = h5.root.track_coords.read()

plt.figure()
plt.title("x vs. px")
plt.plot(tracks[0,0,:], tracks[0,1,:], '.', label="particle 1", ms=1)
plt.plot(tracks[2,0,:], tracks[2,1,:], '.', label="particle 2", ms=1)
plt.plot(tracks[4,0,:], tracks[4,1,:], '.', label="particle 3", ms=1)
plt.xlabel("x")
plt.ylabel("px")
plt.tight_layout()
plt.legend(loc='best')

le01 = le(tracks[0, 0:6, :].transpose(), tracks[1, 0:6, :].transpose())
le23 = le(tracks[2, 0:6, :].transpose(), tracks[3, 0:6, :].transpose())
le45 = le(tracks[4, 0:6, :].transpose(), tracks[5, 0:6, :].transpose())

plt.figure()
plt.title(r'$\ln \frac{\| w(t) \|}{\| w(0) \|}$')
plt.plot(le01, label='LE particle 1')
plt.plot(le23, label='LE particle 2')
plt.plot(le45, label='LE particle 3')
plt.xlabel("turn")
plt.ylabel("LE")
plt.legend(loc='best')

plt.show()
