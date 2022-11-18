#!/usr/bin/env python

import sys, os
import numpy as np
import synergia

lattice = synergia.lattice.MadX_reader().get_lattice('iota', "t1_iota_8_4.madx")
refpart = lattice.get_reference_particle()
print('energy: ', refpart.get_total_energy())
print('momentum: ', refpart.get_momentum())
print('gamma: ', refpart.get_gamma())
print('beta: ', refpart.get_beta())

# The RF cavities have to be tuned (have frequency set) so that the map
# is generated with the correct longitudinal components
closed_orbit = synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
map = synergia.simulation.Lattice_simulator.get_linear_one_turn_map(lattice)

print('map')
print(np.array2string(map, max_line_width=200))

(tunex, tuney, cdt) = synergia.simulation.Lattice_simulator.calculate_tune_and_cdt(lattice)
print('tunex: ', tunex)
print('tuney: ', tuney)

print()
chromaticities = synergia.simulation.Lattice_simulator.get_chromaticities(lattice)

#print('dir(chromaticities): ', dir(chromaticities))
print('chromaticity (x): ', chromaticities.horizontal_chromaticity)
print('chromaticity (y): ', chromaticities.vertical_chromaticity)
