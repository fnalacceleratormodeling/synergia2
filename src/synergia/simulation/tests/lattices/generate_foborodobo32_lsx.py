#!/usr/bin/env python
import synergia
import numpy as np

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0,0]+csmap[1,1])
    asinmu = 0.5*(csmap[0,0]-csmap[1,1])

    if abs(cosmu) > 1.0:
        raise RuntimeError("map is unstable")

    mu =np.arccos(cosmu)

    # beta is positive
    if csmap[0,1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0,1]/np.sin(mu)
    alpha = asinmu/np.sin(mu)
    tune = mu/(2.0*np.pi)

    return (alpha, beta, tune)

#######################################################

lattice = synergia.lattice.Mad8_reader().get_lattice("model", "foborodobo32.lat")

# set frequency of RF cavities
lattice_length = lattice.get_length()          
reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()

# set rf cavity frequency based on length and harmonic number

# harmonic number of this lattice is 128
harmno = 128

# the stable frequency is harmno * beta * c/ring_length
freq = harmno * beta * synergia.foundation.pconstants.c/lattice_length * 1.0e-6
rfwavelen = lattice_length/harmno

print("rf frequency: ", freq)

synergia.utils.write_lsexpr_file(lattice.as_lsexpr(), "foborodobo32_lattice.lsx")

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, 1)

map = synergia.optics.one_turn_map.linear_one_turn_map(lattice_simulator)

[ax, bx, qx] = map2twiss(map[0:2,0:2])
[ay, by, qy] = map2twiss(map[2:4, 2:4])
[az, bz, qz] = map2twiss(map[4:6,4:6])

print("Lattice parameters (assuming uncoupled map)")
print("alpha_x: ", ax, " alpha_y: ", ay)
print("beta_x: ", bx, " beta_y: ", by)
print("q_x: ", qx, " q_y: ", qy)
print("beta_z: ", bz)

print("one turn map")
print(np.array2string(map, max_line_width=300,precision=16))
