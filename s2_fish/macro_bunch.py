#!/usr/bin/env python

from macro_bunch_store import Macro_bunch_store
import Numeric
from mpi4py import MPI

import os.path
import tables

class Macro_bunch:
    def __init__(self):
        self.complete = 0
        self.store = None
        
    def init_test(self,num_per_side,edge_length=1.0):
        '''Paricles uniformly distributed in a cube of size "size"'''
        size = (edge_length,edge_length,edge_length)
        offset = (0.0,0.0,0.0)
        num = (num_per_side,num_per_side,num_per_side)
        local_num = num[0]*num[1]*num[2]
        total_num = local_num # fix me for mpi
        total_current = 1.0
        self.units = Numeric.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = Numeric.array([0.0,0.0,0.0,0.0,0.0,1.1],'d')
        self.particles = Numeric.zeros((7,local_num),'d')
        is_fixedz = 1
        index = 0
        for i in range(0,num[0]):
            for j in range(0,num[1]):
                for k in range(0,num[2]):
                    self.particles[0,index] = offset[0] - size[0]/2.0 + \
                                             size[0]/num[0]*(i+0.5)
                    self.particles[2,index] = offset[1] - size[1]/2.0 + \
                                             size[1]/num[1]*(j+0.5)
                    self.particles[4,index] =  offset[2] - size[2]/2.0 + \
                                             size[2]/num[2]*(k+0.5)
                    self.particles[6,index] = index+1
                    index += 1 
        self.store = Macro_bunch_store(self.particles,local_num,total_num,
                                       total_current,self.units,
                                       self.ref_particle,is_fixedz)
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
                                       1)
        
    def write_particles(self,filename,compress_level=1):
        h5filename = os.path.splitext(filename)[0] + '.h5'
        if MPI.rank == 0:
            f = tables.openFile(h5filename,mode = "w")
            filter = tables.Filters(complevel=compress_level)
            atom = tables.Atom(dtype='Float64',shape=(7,0))
            earray = f.createEArray(f.root,'particles',atom,'Float64',
                                    filters = filter)
            # THIS IS NOT CORRECT. DO NOT LEAVE IT HERE!!!!!
            particles = self.particles
            earray.append(particles)
            for proc in xrange(1,MPI.size):
                parts = MPI.WORLD.Recv(source=proc)
                if parts.shape[1] > 0:
                    earray.append(parts)
            f.close()
        else:
            MPI.WORLD.Send(self.store.get_local_particles(),dest=0)
