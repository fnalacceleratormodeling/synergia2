#!/usr/bin/env python3
import sys, os
import numpy as np
import synergia

def main():

    momentum = 500.0
    four_momentum = synergia.foundation.Four_momentum(synergia.foundation.pconstants.mp)
    four_momentum.set_momentum(momentum)
    refpart = synergia.foundation.Reference_particle(1, four_momentum)

    macroparticles = 100
    realparticles = 5e10

    commxx = synergia.utils.Commxx()

    lattice = synergia.lattice.Lattice("foo")
    lattice.append(synergia.lattice.Lattice_element("marker", "m0"))
    lattice.set_reference_particle(refpart)

    sim = synergia.simulation.Bunch_simulator.create_bunch_train_simulator(
        refpart, 1000, 0.5e10, 1, 1.0)

    sim.set_longitudinal_boundary(
        synergia.bunch.LongitudinalBoundary.
        periodic, 1.0)

    # is it really set?
    bunch = sim.get_bunch(0, 0)
    bb = bunch.get_longitudinal_boundary()
    print('bunch boundary: ', bb[0],': ', bb[1])

    bunch.checkout_particles()
    particles = bunch.get_particles_numpy()
    localnum = bunch.get_local_num()
    print('particles.shape: ', particles.shape)
    print('num local particles: ', localnum)
    z0 = -1
    z1 = +1
    L = (z1 - z0)
    print('Length: ', L)
    dz = L/(localnum-1)
    z = -1
    for i in range(localnum):
        particles[i, 4] = z
        z = z + dz
    bunch.checkin_particles()

    sim.reg_diag_per_turn(synergia.bunch.Diagnostics_particles("foo.h5"), bunch_idx=0)

    sim.set_max_turns(1)

    sc_ops = synergia.collective.Dummy_CO_options()
    stepper = synergia.simulation.Split_operator_stepper(sc_ops, 1)

    print('lattice length: ', lattice.get_length())
    propagator = synergia.simulation.Propagator(lattice, stepper)

    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_TURN)

    propagator.propagate(sim, simlog, 1)

if __name__ == "__main__":
    main()
