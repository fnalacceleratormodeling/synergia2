#!/usr/bin/env python

import synergia

class Pickle_helper:
    __getstate_manages_dict__ = 1
    def __init__(self, *args):
        self.args = args
    def __getinitargs__(self):
        return self.args
    def __getstate__(self):
        return self.__dict__
    def __setstate__(self, state):
        self.__dict__ = state

class Ramp_actions(synergia.simulation.Propagate_actions, Pickle_helper):
    def __init__(self, last_energy):
        synergia.simulation.Propagate_actions.__init__(self)
        Pickle_helper.__init__(self, last_energy)
        self.last_energy = last_energy

    def turn_end_action(self, stepper, bunch, turn_num):
        refpart_energy = bunch.get_reference_particle().get_total_energy()
        stepper.get_lattice_simulator().get_lattice().get_reference_particle().set_total_energy(refpart_energy)
        self.last_energy = refpart_energy
        stepper.get_lattice_simulator().update()
