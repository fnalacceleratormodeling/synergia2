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
    def __init__(self, sc):
        synergia.simulation.Propagate_actions.__init__(self)
        Pickle_helper.__init__(self, sc)
        self.sc = sc
    def turn_end_action(self, stepper, bunch, turn_num):
        print "telling sc to dump", "ex_%04d.h5" %turn_num
        self.sc.set_files("ex_%04d.h5" %turn_num, "ey_%04d.h5" %turn_num)

