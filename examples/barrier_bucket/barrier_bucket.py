#!/usr/bin/python

import sys, os
import numpy as np
import synergia
import matplotlib.pyplot as plt
import h5py

if len(sys.argv) != 2:
    print("usage: barrier_bucket.py <num-steps>")
    sys.exit(10)

num_steps = int(sys.argv[1])
print("using ", num_steps, "steps")

lattice = synergia.lattice.Lattice("foo", synergia.lattice.MadX_adaptor_map())

drlen = 10.0
dr = synergia.lattice.Lattice_element("drift", "dr")
dr.set_double_attribute("l", drlen)

lattice.append(dr)

# # betagamma=3/4, gamma = 5/4, beta = 3/5
# beta = 0.6

mp = synergia.foundation.pconstants.mp
energy = 0.8 + mp
refpart = synergia.foundation.Reference_particle(1, mp, energy)
four_momentum = refpart.get_four_momentum()

lattice.set_reference_particle(refpart)

print(lattice.as_string())

print("particle created")
print('Energy: ', refpart.get_total_energy())
print('Momentum: ', refpart.get_momentum())
print('gamma: ', refpart.get_gamma())
print('beta: ', refpart.get_beta())

beta_2 = refpart.get_beta() * 10.0/9.0
four_momentum_2 = synergia.foundation.Four_momentum(mp)
four_momentum_2.set_beta(beta_2)

#dpop = four_momentum_2.get_momentum()/four_momentum.get_momentum() - 1.0
dpop = 0.1
print('dpop: ', dpop)

commxx = synergia.utils.Commxx()
bunch = synergia.bunch.Bunch(refpart, 2, 1.0e9, commxx)
bunch.set_bucket_index(0)
bunch.set_bucket_barrier_length(drlen)

local_particles = bunch.get_local_particles()
local_particles[0, 5] = dpop

bunch_simulator = synergia.simulation.Bunch_simulator(bunch)
bunch_simulator.add_per_step(synergia.bunch.Diagnostics_bulk_track("pystep_tracks.h5", 2))
bunch_simulator.add_per_turn(synergia.bunch.Diagnostics_bulk_track("pyturn_tracks.h5", 2))

stepper = synergia.simulation.Independent_stepper(lattice, 1, num_steps)
propagator = synergia.simulation.Propagator(stepper)

propagator.propagate(bunch_simulator, 200, 200, 1)

del propagator
del stepper
del bunch_simulator
del bunch

sys.exit(0)

h5 = h5py.File('tracks.h5', 'r')
s = h5.get('s')
track_coords = h5.get('track_coords')
# track_coords is <#turns> x <#particles> x <coordinate>
ax1 = plt.subplot(211)

ax1.plot(s, track_coords[:, 0, 4], label='cdt')
xticklabels1 = ax1.get_xticklabels()
plt.setp(xticklabels1, visible=False)
plt.ylabel('cdt position')
plt.legend(loc='best')

ax2 = plt.subplot(212, sharex=ax1)
ax2.plot(s, track_coords[:, 0, 5], label='dp/p')
plt.xlabel('s position [m]')
plt.ylabel('dp/p')
plt.legend(loc='best')

plt.subplots_adjust(hspace=.001)

plt.show()

        
