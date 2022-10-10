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
    # The arguments to __init__ are what the Ramp_actions instance is
    # initialized with
    def __init__(self, numbers):
        synergia.simulation.Propagate_actions.__init__(self)
        # pickling the arguments to the initializer allows the
        # module to resume after checkpointing.  They should be in the same
        # order as the arguments to __init__.
        Pickle_helper.__init__(self, numbers)
        self.numbers = numbers

    def turn_end_action(self, stepper, bunch, turn_num):
        if bunch.get_comm().get_rank() == 0:
            print('turn number: ', turn_num, ', number: ', self.numbers[turn_num],', reference particle state: ', bunch.get_reference_particle().get_state())
            print('lattice momentum: ', stepper.get_lattice_simulator().get_lattice().get_reference_particle().get_momentum())
            print('bunch momentum: ', bunch.get_reference_particle().get_momentum())
            print('bunch design momentum: ', bunch.get_design_reference_particle().get_momentum())
            stepper.get_lattice_simulator().get_lattice().get_reference_particle().get_four_momentum().set_momentum(
                bunch.get_reference_particle().get_four_momentum().get_momentum())
            bunch.get_design_reference_particle().get_four_momentum().set_momentum(
                bunch.get_reference_particle().get_four_momentum().get_momentum())
            print()
        stepper.get_lattice_simulator().get_lattice().get_reference_particle().get_four_momentum().set_momentum(bunch.get_reference_particle().get_four_momentum().get_momentum())
            
    # other possible actions which could be present.
    # actions that are not present will default to internal stubs
    def first_action(self, stepper, bunch):
        pass

    def step_end_action(self, stepper, step, bunch, turn_num, step_num):
        pass

    def before_resume_action(self, stepper, bunch):
        pass
