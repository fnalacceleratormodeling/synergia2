#!/usr/bin/env python

import pytest
import numpy
import synergia

@pytest.fixture
def prop_fixture():
    fodo_madx = """
beam, particle=proton,pc=3.0;

o: drift, l=8.0;
f: quadrupole, l=2.0, k1=0.071428571428571425, a1=1.0;
d: quadrupole, l=2.0, k1=-0.071428571428571425;

fodo: sequence, l=20.0, refer=entry;
fodo_1: f, at=0.0;
fodo_2: o, at=2.0;
fodo_3: d, at=10.0;
fodo_4: o, at=12.0;
endsequence;
"""
    reader = synergia.lattice.MadX_reader()
    reader.parse(fodo_madx)
    lattice = reader.get_lattice('fodo')
    lattice.set_all_string_attribute('extractor_type', 'libff')
    stepper = synergia.simulation.Independent_stepper_elements(1)
    propagator = synergia.simulation.Propagator(lattice, stepper)

    return {'lattice': lattice, 'propagator': propagator}

def test_create_prop(prop_fixture):
    return True

def test_prop_get_lattice_elements(prop_fixture):
    nlattice_elem = len(prop_fixture['lattice'].get_elements())
    nprop_latt_elem = len(prop_fixture['propagator'].get_lattice_elements())
    assert nlattice_elem == nprop_latt_elem
    return True

def test_prop_get_lattice(prop_fixture):
    nlattice_elem = len(prop_fixture['lattice'].get_elements())
    nprop_latt_elem = len(prop_fixture['propagator'].get_lattice().get_elements())
    assert nlattice_elem == nprop_latt_elem
    return True

def test_modify_lattice_elements(prop_fixture):
    # print(prop_fixture['propagator'].get_lattice().get_elements()[0])
    orig_a1 = prop_fixture['propagator'].get_lattice().get_elements()[0].get_double_attribute('a1')
    new_a1 = orig_a1 + 100.0
    prop_fixture['propagator'].get_lattice().get_elements()[0].set_double_attribute('a1', new_a1)
    assert prop_fixture['propagator'].get_lattice().get_elements()[0].get_double_attribute('a1') == new_a1

    slices = list(prop_fixture['propagator'].get_lattice_element_slices())
    assert slices[0].get_lattice_element().get_double_attribute('a1') == new_a1
    return 0

def test_modify_reference_particle_energy(prop_fixture):
    orig_energy = prop_fixture['propagator'].get_lattice().get_reference_particle().get_total_energy()
    new_energy = orig_energy + 2.0
    prop_fixture['propagator'].get_lattice().get_reference_particle().set_total_energy(new_energy)
    assert prop_fixture['propagator'].get_lattice().get_reference_particle().get_total_energy() == new_energy
    return True

def test_set_get_checkpoint(prop_fixture):
    init_period = prop_fixture['propagator'].get_checkpoint_period()
    cp_period = init_period + 1
    prop_fixture['propagator'].set_checkpoint_period(cp_period)
    assert cp_period == prop_fixture['propagator'].get_checkpoint_period()
    return True

def test_set_get_final_checkpoint(prop_fixture):
    init_fcp = prop_fixture['propagator'].get_final_checkpoint()
    prop_fixture['propagator'].set_final_checkpoint(not init_fcp)
    assert prop_fixture['propagator'].get_final_checkpoint() is not init_fcp
    return True
