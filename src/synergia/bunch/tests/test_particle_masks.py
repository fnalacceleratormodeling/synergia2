#!/usr/bin/env python

import pytest
import numpy as np
import synergia

mass = 1.0
betagamma = 3/4
gamma = 5/4
energy = mass*gamma

macroparticles=16
realparticles=5e10


#@pytest.fixture
def bunch_fixture():
    refpart = synergia.foundation.Reference_particle(1, mass, energy)
    bunch = synergia.bunch.Bunch(refpart, macroparticles, realparticles, synergia.utils.Commxx())
    bunch.checkin_particles()
    return bunch

def test_mask_access(bunch_fixture):
    bunch_fixture.checkout_particles()
    masks = bunch_fixture.get_particle_masks_numpy()

def test_mask_initial_state(bunch_fixture):
    # all masks should be 1 when bunch is created
    bunch_fixture.checkout_particles()
    masks = bunch_fixture.get_particle_masks_numpy()
    assert masks.sum() == macroparticles

# for running from the command line, comment out pytest.fixture
if __name__ == "__main__":
    bf = bunch_fixture()

    test_mask_access(bf)
    test_mask_initial_state(bf)
