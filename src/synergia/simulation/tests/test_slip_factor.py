#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles=16
realparticles=4.0e10
nturns=1
mp = synergia.foundation.pconstants.mp
dpp = 1.0e-3

# prop_fixture is a propagator
#@pytest.fixture
def prop_fixture():

    channel_madx0 = """
beam, particle=proton,energy=pmass+0.8;
b: sbend, angle=pi/2, l=pi/2;
channel: sequence, l=pi/2, refer=entry;
b, at=0.0;
endsequence;
"""
    channel_madx1 = """
beam, particle=proton,energy=pmass+0.8;
d: drift, l=1.0;
b: sbend, angle=pi/2, l=pi/2;
channel: sequence, l=4.0+2*pi, refer=entry;
d, at=0.0;
b, at=1.0;
d, at=1.0+pi/2;
b, at=2.0+pi/2;
d, at=2.0+pi;
b, at=3.0+pi;
d, at=3.0+3*pi/2;
b, at=4.0+3*pi/2;
endsequence;
"""
    channel_mad2 = """
beam, particle=proton,energy=pmass+0.8;
d: drift, l=1.0;
b: sbend, angle=pi/2, l=pi/2;channel: line=(b, d, b, d, b, d, b, d);
"""
    channel_madx = """
beam, particle=proton,energy=pmass+0.8;
d: drift, l=1.0;
b: sbend, angle=pi/2, l=pi/2;
channel: sequence, l=4.0+2*pi, refer=entry;
b, at=0.0;
d, at=pi/2;
b, at=1.0+pi/2;
d, at=1.0+pi;
b, at=2.0+pi;
d, at=3.0+3*pi/2;
b, at=3.0+3*3*pi/2;
d, at=3.0+2*pi;
endsequence;
"""
    channel_madx5 = """
beam, particle=proton,energy=pmass+0.8;
d: drift, l=1.0;
b: sbend, angle=pi/2, l=pi/2;
channel: sequence, l=pi/2, refer=entry;
b, at=0.0;
endsequence;
"""


    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice('channel')
    lattice.set_all_string_attribute('extractor_type', 'libff')
    synergia.simulation.Lattice_simulator.tune_linear_lattice(lattice)
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator

# The beamline consists of a four pi/2 bend. Two particles are propagated, one at the design momentum
# and one at an offset momemtum. The particles propagate in a circular orbit but with different radii. We
# can calculate the difference in path length and transit time.

# p1: design momentum
# p2: offset momentum = (1+delta)*p1
# R1: radius of particle p1 = B/(e p1)
# R2: radius of particle p2 = B/(e p2) = B/(e p1 (1+dpp))
#
# cos (particle2 angle) = 1 - R1/R2
#                       = delta/(1 + delta)
#
# Design path length: pi/2
# offset momentum path length = arccos(delta/(1+delta))
# delta path length = arccos(delta/(1+delta)) - pi/2
# c*reftime = pi/2 * (1/beta)


def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(ref_part, macroparticles, realparticles)
    bunch = sim.get_bunch()
    bunch.checkout_particles()

    energy = ref_part.get_total_energy()

    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
 
    bunch.checkin_particles()
    return sim


def test_slip_factor(prop_fixture):

    #print('lattice: ')
    #print(prop_fixture.get_lattice())

    refpart = prop_fixture.get_lattice().get_reference_particle()

    #synergia.simulation.Lattice_simulator.set_closed_orbit_tolerance(1.0/524288)
    #assert synergia.simulation.Lattice_simulator.get_closed_orbit_tolerance() == 1.0/524288

    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())
    bunch = sim.get_bunch(0, 0)
    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False)
    prop_fixture.propagate(sim, simlog, nturns)
 
    beta1 = bunch.get_reference_particle().get_beta()
    gamma1 = bunch.get_reference_particle().get_gamma()
    betagamma1 = beta1*gamma1

    betagamma2 = (1+dpp) * betagamma1
    gamma2 = np.sqrt(1+betagamma2**2)
    beta2 = betagamma2/gamma2

    lattice = prop_fixture.get_lattice()
    bend = None
    drift = None
    for elem in lattice.get_elements():
        #print(elem)
        if not bend and elem.get_type() == synergia.lattice.element_type.sbend:
            #print('found bend')
            bend = elem
            continue
        elif not drift and elem.get_type() == synergia.lattice.element_type.drift:
            #print('found drift')
            drift = elem
            continue
        if bend and drift:
            break

    if not bend:
        raise RuntimeError('no bend element found')
    
    # L = R theta => R = L/theta
    # R is radius of curvature of bend magnet ~ p/eB
    # R2/R1 = p2/p1
    R1 = bend.get_length()/bend.get_bend_angle()
    R2 = R1 * (1+dpp)

    L1 = R1 * np.pi/2
    L2 = R2 * np.arccos(1 - R1/R2)

    if drift:
        Ld = drift.get_length()
    else:
        Ld = 0.0

    # compaction factor [(L2 - L1)/L1] / dpp
    # The ring is designed so that the orbit in the drifts
    # is parallel to the axis so the differing momentum only
    # changes path length in the bends*
    # with momentum
    compact_calc = (L2/L1 - 1)/dpp

    cT1 = L1/beta1
    cT2 = L2/beta2
    cTdrift1 = Ld/beta1
    cTdrift2 = Ld/beta2

    # slip factor = [(T2 - T1)/T1] / dpp
    slip_calc = ( (cT2 - cT1 + cTdrift2 - cTdrift1)/(cT1 + cTdrift1) )/ dpp

    # The drifts should have the same length for both
    # momentum particles since the orbit is parallel to
    # the axis.
    print('compaction_calc', compact_calc)
    print('slip_calc: ', slip_calc)
    print('slip_factor + 1/gamma**2: ', slip_calc + 1/gamma1**2)

    bunch = sim.get_bunch(0, 0)
    bunch.checkout_particles()
    part = bunch.get_particles_numpy()
    part[1, 0] = (R2-R1)
    part[1, 5] = dpp
    bunch.checkin_particles()

    print('initial particles')
    print(part[0,:])
    print(part[1,:])
    print()
    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False)
    prop_fixture.propagate(sim, simlog, nturns)

    print('final particles')
    bunch.checkout_particles()
    print(part[0, :])
    print(part[1, :])

    lattice.get_reference_particle().set_state(dpp, 0.0, 0.0, 0.0, 0.0, dpp)

    synergia.simulation.Lattice_simulator.set_closed_orbit_tolerance(1.0/65536)
    assert synergia.simulation.Lattice_simulator.get_closed_orbit_tolerance() == 1.0/65536
    chrom = synergia.simulation.Lattice_simulator.get_slip_factors(lattice)
    slip_factor = chrom.slip_factor
    momentum_compact = chrom.momentum_compaction
    print('lattice_simulator slip factor: ', slip_factor)
    print('lattice simulator momentum compaction: ', momentum_compact)

def main():
    pf = prop_fixture()
    test_slip_factor(pf)

if __name__ == "__main__":
    main()
