#!/usr/bin/env python
import sys, os
import numpy as np
import synergia
import pytest

macroparticles=16
realparticles=4.0e10
# lag 1/120 is a phase angle of 2pi/120 or pi/60 or 3 degrees
# V = 0.2 MV * sin(pi/60) = 
expected_delta_E = 0.0002*np.sin(np.pi/60)
print('expected delta E/turn: ', expected_delta_E)
nturns=1
dpop_offset = 1.0e-3
transmom_offset = 0.001
mp = synergia.foundation.pconstants.mp

# prop_fixture is a propagator
@pytest.fixture
def prop_fixture():
    
    # rfcavity not accelerating
    channel_madx = """
beam, particle=proton,pc=0.75*pmass;
rfc: rfcavity, l=0.0, volt=0.0, harmon=1;
q: quadrupole, l=1, k1=0.0625;
channel: sequence, l=20.0, refer=centre;
rfc,volt=0, at=20.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(channel_madx)
    lattice = reader.get_lattice('channel')
    lattice.set_all_string_attribute('extractor_type', 'libff')
    synergia.simulation.Lattice_simulator.tune_circular_lattice(lattice)
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return propagator



def create_simulator(ref_part):
    sim = synergia.simulation.Bunch_simulator.create_single_bunch_simulator(ref_part, macroparticles, realparticles)
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    lp[:, 0:6] = 0.0
    bunch.checkin_particles()
    return sim


def test_accel1(prop_fixture):

    print('lattice: ')
    print(prop_fixture.get_lattice())

    lattice = prop_fixture.get_lattice()
    frequency = synergia.simulation.Lattice_simulator.get_rf_frequency(lattice)

    # set up particle to get momentum increase as if it had been 
    # accelerated
    ref_part = lattice.get_reference_particle()
    e = ref_part.get_total_energy()
    p = ref_part.get_momentum()
    enew = e + expected_delta_E
    pnew = np.sqrt(enew**2 - mp**2)
    expected_dpop = (pnew/p)-1
    # for expected kick, want V sin(w/c * cdt) = expected_de
    # so cdt = arcsin( expected_de/v) * c/w
    new_cdt = (np.pi/60) * synergia.foundation.pconstants.c/frequency
    print('expected dp/p: ', expected_dpop)
    print('cdt for expected dp/p: ', new_cdt)
    sim = create_simulator(prop_fixture.get_lattice().get_reference_particle())
    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()

    lp[0, 4] = new_cdt
    bunch.checkin_particles()

    #simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.INFO_TURN, False)
    simlog = synergia.utils.parallel_utils.Logger(0, synergia.utils.parallel_utils.LoggerV.DINFO, True)
    print('before propagation, lp[0,4]: ', lp[0,4])
    print('before propagation ', lp[0,5])
    prop_fixture.propagate(sim, simlog, nturns)

    bunch = sim.get_bunch()
    bunch.checkout_particles()
    lp = bunch.get_particles_numpy()
    assert lp[0, 5] == expected_dpop

    # cavity at the end of the line, so cdt should remain 0
    assert lp[0,4] == pytest.approx(0.0)
    print('lp[0]: ', lp[0,:])

    #assert False

    

def main():
    pf = prop_fixture()
    test_accel1(pf)

if __name__ == "__main__":
    main()
