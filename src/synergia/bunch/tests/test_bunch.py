#!/usr/bin/env python

import pytest
import numpy
import synergia

mass = 1.0
betagamma = 3/4
gamma = 5/4
energy = mass*gamma

@pytest.fixture
def bunch_fixture():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    macroparticles=10
    realparticles=5e10
    bunch = synergia.bunch.Bunch(refpart, macroparticles, realparticles, synergia.utils.Commxx())
    return bunch

def test_bunch_1(bunch_fixture):
    assert bunch_fixture.get_total_num() == 10
    assert bunch_fixture.get_real_num() == pytest.approx(5e10)
    bunch_fixture.set_real_num(6e10)
    assert bunch_fixture.get_real_num() == pytest.approx(6e10)


