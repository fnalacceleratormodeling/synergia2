#!/usr/bin/env python

from macro_bunch_store import Macro_bunch_store
import Numeric

class Macro_bunch:
    def __init__(self):
        self.complete = 0

    def init_test(self):
        size = (1.0,1.0,1.0)
        offset = (0.0,0.0,0.0)
        num = (10,10,10)
        local_num = num[0]*num[1]*num[2]
        total_num = local_num # fix me for mpi
        total_current = 1.0
        self.units = Numeric.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = Numeric.array([0.0,0.0,0.0,0.0,0.0,1.1],'d')
        self.particles = Numeric.array((7,local_num),'d')
        is_z = 0
        id = 0
        for i in range(0,num[0]):
            for j in range(0,num[1]):
                for k in range(0,num[2]):
                    id += 1
                    self.particles[(0,id)] = offset[0] - size[0]/2.0 + \
                                             size[0]/(num[0]-1)*(i+0.5)
                    self.particles[(2,id)] = offset[1] - size[1]/2.0 + \
                                             size[1]/(num[1]-1)*(j+0.5)
                    self.particles[(4,id)] =  offset[2] - size[2]/2.0 + \
                                             size[2]/(num[2]-1)*(k+0.5)
                    self.particles[(6,id)] = id
        self.store = Macro_bunch_store(self.particles,local_num,total_num,
                                       total_current,self.units,
                                       self.ref_particle,is_z)
        self.complete = 1
        
    def init_from_beambunch(self, beambunch):
        pass
    
        
