#!/usr/bin/env python

import synergia
import memory

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

class Trackmem_actions(synergia.simulation.Propagate_actions, Pickle_helper):
    def __init__(self):
        synergia.simulation.Propagate_actions.__init__(self)
        Pickle_helper.__init__(self)
        self.file = open('trackmem.log','w')
    def turn_end_action(self, stepper, bunch, turn_num):
        mem = memory.memory()
        res = memory.resident()
        stack = memory.stacksize()
        print "trackmem turn %d: %d %d %d" %(turn_num,mem,res,stack)
        self.file.write("%d %d %d %d\n" %(turn_num,mem,res,stack))
