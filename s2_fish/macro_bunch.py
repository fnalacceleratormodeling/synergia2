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
        self.particles = Numeric.zeros((7,local_num),'d')
        is_z = 0
        index = 0
        for i in range(0,num[0]):
            for j in range(0,num[1]):
                for k in range(0,num[2]):
                    self.particles[0,index] = offset[0] - size[0]/2.0 + \
                                             size[0]/(num[0]-1)*(i+0.5)
                    self.particles[2,index] = offset[1] - size[1]/2.0 + \
                                             size[1]/(num[1]-1)*(j+0.5)
                    self.particles[4,index] =  offset[2] - size[2]/2.0 + \
                                             size[2]/(num[2]-1)*(k+0.5)
                    self.particles[6,index] = index+1
                    index += 1 
        self.store = Macro_bunch_store(self.particles,local_num,total_num,
                                       total_current,self.units,
                                       self.ref_particle,is_z)
        self.complete = 1
        
    def init_from_bunch(self, bunch):
        (Cxy, Cxpyp, Cz, Czp) = bunch.beam_parameters.get_conversions()
        self.units = Numeric.array([Cxy,Cxpyp,Cxy,Cxpyp,Cz,Czp],'d')
        self.store = Macro_bunch_store(bunch.particles(),
                                       bunch.num_particles_local(),
                                       bunch.num_particles(),
                                       bunch.current(),
                                       self.units,
                                       bunch.reference_particle(),
                                       0)
    
        
