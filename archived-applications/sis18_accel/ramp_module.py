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

    # def step_end_action(self, stepper, step, bunch, turn_num, step_num):
    #     bunch_energy = bunch.get_reference_particle().get_total_energy()
    #     if bunch_energy != self.last_energy:
    #         print "step ", step_num, " energy changed from ", self.last_energy, " to ", bunch_energy
    #         stepper.get_lattice_simulator().get_lattice().get_reference_particle().set_total_energy(bunch_energy)
    #         self.last_energy = bunch_energy
    #         stepper.get_lattice_simulator().update()
    #         sliced_beamline = stepper.get_lattice_simulator().get_chef_lattice().get_sliced_beamline()
    #         print "sliced_beamline energy: ", sliced_beamline.Energy()
