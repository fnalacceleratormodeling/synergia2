#!/usr/bin/env python

from macro_bunch_store import Macro_bunch_store
import Numeric
import RandomArray
from mpi4py import MPI

import os.path
import tables
import time
import math

class Macro_bunch:
    def __init__(self):
        self.complete = 0
        self.particles = None
        self.local_num = None
        self.total_num = None
        self.total_current = None
        self.is_fixedz = None
        self.ref_particle = None

    def get_local_particles(self):
        return self.particles
    
    def get_num_particles_local(self):
        return self.local_num
    
    def get_store(self):
        if self.particles:
            return Macro_bunch_store(self.particles,
                                     self.local_num,
                                     self.total_num,
                                     self.total_current,
                                     self.units,
                                     self.ref_particle,
                                     self.is_fixedz)
        else:
            return None
        
    def convert_to_fixedz(self):
        self.get_store().convert_to_fixedz()
        self.is_fixedz = 1

    def convert_to_fixedt(self):
        self.get_store().convert_to_fixedt()
        self.is_fixedz = 0
        
    def init_test(self,num_per_side,edge_length=1.0):
        '''Particles uniformly distributed in a cube of size "size"'''
        size = (edge_length,edge_length,edge_length)
        offset = (0.0,0.0,0.0)
        num = (num_per_side,num_per_side,num_per_side)
        local_num = num[0]*num[1]*num[2]
        total_num = local_num # fix me for mpi
        total_current = 1.0
        self.units = Numeric.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = Numeric.array([0.0,0.0,0.0,0.0,0.0,1.1],'d')
        self.particles = Numeric.zeros((7,local_num),'d')
        self.is_fixedz = 1
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
        self.local_num = local_num
        self.total_num = total_num
        self.total_current = total_current
        self.complete = 1

    def init_sphere(self,num,radius):
        '''Particles uniformly distributed in a sphere of radius "radius"'''
        RandomArray.seed(17+MPI.rank*51,59+MPI.rank*23)
        offset = (0.0,0.0,0.0)
        local_num = num/MPI.size #jfa: could be more precise...
        total_num = MPI.WORLD.Allreduce(local_num,MPI.SUM)
        total_current = 1.0
        self.units = Numeric.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = Numeric.array([0.0,0.0,0.0,0.0,0.0,1.1],'d')
        self.particles = Numeric.zeros((7,local_num),'d')
        self.is_fixedz = 1
        index = 0
        added = 0
        discarded = 0
        chunk_size = 1000
        t0 = time.time()
        p = (RandomArray.random([6,chunk_size])-0.5)*2.0*radius
        index = 0
        while added < local_num:
            if ((p[0,index]**2 + p[2,index]**2 + p[4,index]**2) < radius**2):
                self.particles[0:6,added] = p[:,index]
                self.particles[6,added] = added + 1
                added += 1
            else:
                discarded += 1
            index += 1
            if index >= chunk_size:
                p = (RandomArray.random([6,chunk_size])-0.5)*2.0*radius
                index = 0
        t1 = time.time()
#        print "pi =",6.0*added/(1.0*added+discarded),"in",t1-t0,"secs"
        self.local_num = local_num
        self.total_num = total_num
        self.total_current = total_current
        self.complete = 1

    def init_cylinder(self,num,radius,length):
        '''Particles uniformly distributed in a cylinder of radius "radius"
        and length 2pi'''
        offset = (0.0,0.0,0.0)
        local_num = num
        total_num = local_num # fix me for mpi
        total_current = 1.0
        self.units = Numeric.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = Numeric.array([0.0,0.0,0.0,0.0,0.0,1.1],'d')
        self.particles = Numeric.zeros((7,local_num),'d')
        self.is_fixedz = 1
        index = 0
        added = 0
        discarded = 0
        chunk_size = 1000
        t0 = time.time()
        p = (RandomArray.random([6,chunk_size])-0.5)*2.0*radius
        index = 0
        while added < local_num:
            if ((p[0,index]**2 + p[2,index]**2) < radius**2):
                self.particles[0:6,added] = p[:,index]
                self.particles[4,added] *= length/(2.0*radius)
                self.particles[6,added] = added + 1
                added += 1
            else:
                discarded += 1
            index += 1
            if index >= chunk_size:
                p = (RandomArray.random([6,chunk_size])-0.5)*2.0*radius
                index = 0
        t1 = time.time()
        self.local_num = local_num
        self.total_num = total_num
        self.total_current = total_current
        self.complete = 1
        
    def init_from_bunch(self, bunch):
        (Cxy, Cxpyp, Cz, Czp) = bunch.beam_parameters.get_conversions()
        self.units = Numeric.array([Cxy,Cxpyp,Cxy,Cxpyp,Cz,Czp],'d')
        self.particles = bunch.particles()
        self.local_num = bunch.num_particles_local()
        self.total_num = bunch.num_particles()
        self.total_current = bunch.current()
        print "current from bunch =",bunch.current()
        self.ref_particle = bunch.reference_particle()
        self.is_fixedz = 1
        
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
            MPI.WORLD.Send(self.particles,dest=0)

    def write_particles_text(self,filename):
        if MPI.rank == 0:
            f = open(filename,"w")
            for proc in xrange(1,MPI.size):
                parts = MPI.WORLD.Recv(source=proc)
                for i in range(0,parts.shape[1]):
                    f.write("%g %g %g %g %g %g %g\n" % \
                            tuple(parts[:,i]))
            parts = self.particles
            for i in range(0,parts.shape[1]):
                f.write("%g %g %g %g %g %g %g\n" % \
                        tuple(parts[:,i]))
            f.close()
        else:
            MPI.WORLD.Send(self.particles(),dest=0)
