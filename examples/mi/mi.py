#!/usr/bin/env python
import synergia
import numpy as np
from mi_options import opts
from synergia.optics.one_turn_map import linear_one_turn_map

# quick and dirty twiss parameter calculator from 2x2 courant-snyder map array
def map2twiss(csmap):
    cosmu = 0.5 * (csmap[0, 0] + csmap[1, 1])
    asinmu = 0.5 * (csmap[0, 0] - csmap[1, 1])

    if abs(cosmu) > 1.0:
        print "error, map is unstable"
    mu = np.arccos(cosmu)

    # beta is positive
    if csmap[0, 1] < 0.0:
        mu = 2.0 * np.pi - mu

    beta = csmap[0, 1] / np.sin(mu)
    alpha = asinmu / np.sin(mu)
    tune = mu / (2.0 * np.pi)

    return (alpha, beta, tune)

num_macro_particles = opts.num_macro_particles
seed = opts.seed
#grid = [16, 16, 128]
grid = [16, 16, 16]
num_real_particles = 1e11
num_steps = opts.num_steps
num_turns = opts.num_turns
map_order = opts.map_order
norm_emit = opts.norm_emit
stdz = 0.3
harmno = 588

lattice = synergia.lattice.Mad8_reader().get_lattice("ring_p_q605", "mi20-egs-thinrf.lat")
lattice_length = lattice.get_length()
print "lattice length: ", lattice_length
reference_particle = lattice.get_reference_particle()
energy = reference_particle.get_total_energy()
beta = reference_particle.get_beta()
gamma = reference_particle.get_gamma()

print "particle energy: ", energy
print "particle beta: ", beta
print "particle gamma: ", gamma

# set rf cavity frequency
# harmno * beta * c/ring_length
freq = harmno * beta * synergia.foundation.pconstants.c / lattice_length
print "RF freq: ", freq

# rf cavity voltage, is 1.0 MV total distributed over 18 cavities.  MAD8
# expects cavities voltages in  units of MV.
# At injection, the Main Injector is below transition so the phase is 0
for elem in lattice.get_elements():
    if elem.get_type() == "rfcavity":
        elem.set_double_attribute("volt", opts.rf_voltage)
        elem.set_double_attribute("freq", freq)

lattice_simulator = synergia.simulation.Lattice_simulator(lattice, map_order)

map = linear_one_turn_map(lattice_simulator)
print "one turn map from synergia2.5 infrastructure"
print np.array2string(map, max_line_width=200)

[[ax, ay], [bx, by]] = synergia.optics.get_alpha_beta(map)
print "Lattice functions assuming uncoupled map:"
print "alpha x: ", ax
print "alpha y: ", ay
print "beta x: ", bx
print "beta y: ", by

[az, bz, qz] = map2twiss(map[4:6, 4:6])
print "alpha z (better be small): ", az
print "beta z: ", bz

# dummy collective operator
no_op = synergia.simulation.Dummy_collective_operator("stub")

# space charge collective operator
#space_charge = synergia.simulation.Collective_operator("space charge")

#define stepper with space charge (what is the new way?)
#stepper_sc = synergia.simulation.Split_operator_stepper(lattice_simulator, space_charge,
#                                          num_steps)


# define stepper with dummy collective operator
stepper_noop = synergia.simulation.Split_operator_stepper(lattice_simulator, no_op ,
                                          num_steps)

emit = norm_emit / (beta * gamma)
print "generating particles with transverse emittance: ", emit
print "expected std_x: ", np.sqrt(emit * bx)
print "expected std_y: ", np.sqrt(emit * by)
print "expected std_z: ", opts.stdz
print "expected std_dpop: ", opts.stdz / bz

covar = synergia.optics.matching._get_correlation_matrix(map, np.sqrt(emit * bx), np.sqrt(emit * by), opts.stdz, bz)
print "covariance matrix"
print np.array2string(covar, max_line_width=200)

bunch = synergia.optics.generate_matched_bunch(lattice_simulator, np.sqrt(emit * bx), np.sqrt(emit * by), stdz, num_real_particles, num_macro_particles,
                                        seed=seed)



# get particles for inspection and manipulations
particles = bunch.get_local_particles()

print "Generated bunch properties:"
print "     x mean: ", particles[:, 0].mean(), " std: ", particles[:, 0].std()
print "     y mean: ", particles[:, 2].mean(), " std: ", particles[:, 2].std()
print "     z mean: ", particles[:, 4].mean(), " std: ", particles[:, 4].std()

particles[:, 0] = particles[:, 0] + opts.x_offset
particles[:, 2] = particles[:, 2] + opts.y_offset
particles[:, 4] = particles[:, 4] + opts.z_offset

multi_diagnostics_step = synergia.bunch.Multi_diagnostics()
multi_diagnostics_turn = synergia.bunch.Multi_diagnostics()
multi_diagnostics_turn.append(synergia.bunch.Diagnostics_full2(bunch, "mi_full2.h5"))
multi_diagnostics_turn.append(synergia.bunch.Diagnostics_particles(bunch, "mi_particles.h5"))

propagator = synergia.simulation.Propagator(stepper_noop)
propagator.propagate(bunch, num_turns, multi_diagnostics_step,
                      multi_diagnostics_turn, opts.verbose)
