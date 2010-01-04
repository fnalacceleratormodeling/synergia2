#!/usr/bin/env python

import numpy
import os
import tables
import time
import shutil
import operator
from mpi4py import MPI

from job_manager import create_new_directory


class Tracker:
    def __init__(self,root_dir,fractuple,sub_dir=None,dest_dir=".",track_dir='tracks',save_period=20):
        # if fractuple is of the form (a,b), interpret as a/b
        # else interpret as 1/fractuple 
        if operator.isSequenceType(fractuple):
            self.numer = fractuple[0]
            self.denom = fractuple[1]
        else:
            self.numer = 1
            self.denom = fractuple
        if not sub_dir:
            sub_dir = os.environ['USER'] + str(os.getpid())
        if MPI.COMM_WORLD.Get_rank() == 0:
            self.base_dir = create_new_directory(os.path.join(root_dir,
                                                              sub_dir),
                                                 0,overwrite=0)
        else:
            self.base_dir = "/dummy"
        self.dir = os.path.join(self.base_dir,track_dir)
        self.dest_dir = dest_dir
        self.track_dir = track_dir
        self.add_times = []
        self.copy_time = None
        self.open = 1
        self.count = 0
        self.save_period=save_period
        
    def _write_track_text(self,part,s):
        id = int(part[6])
        fname = index2filename(id,self.dir,"dat")
        if not os.path.isdir(os.path.dirname(fname)):
            os.makedirs(os.path.dirname(fname))
        f = open(fname,"a")
        f.write("%g %g %g %g %g %g " % \
                tuple(part[0:6]))
        f.write("%g\n" % s)
        f.close()


    def _write_track(self,part,s):
        id = int(part[6])
        fname = index2filename(id,self.dir,"h5")
        if not os.path.exists(fname):
            if not os.path.isdir(os.path.dirname(fname)):
                os.makedirs(os.path.dirname(fname))
            f = tables.openFile(fname,mode = "w")
            atom = tables.Atom(dtype='Float64',shape=(7,0))
            earray = f.createEArray(f.root,'points',atom,'Float64')
        else:
            f = tables.openFile(fname,mode = "a")
            earray = f.root.points
        point = numpy.zeros((7,1),'d')
        point[0:6,0] = part[0:6]
        point[6,0] = s
        earray.append(point)
        f.close()

    def add(self, bunch, s):
        if not self.open:
            print "Warning: attempt to add to tracker after closing. Nothing added."
            return
        last_flag = -999
        if MPI.COMM_WORLD.Get_rank() == 0:
            t0 = time.time()
            for i in range(0,bunch.get_num_particles_local()):
                id = int(bunch.get_local_particles()[6,i])
                if (self.numer*id % self.denom) < self.numer:
                    self._write_track(bunch.get_local_particles()[:,i],s)            
            for proc in xrange(1,MPI.COMM_WORLD.Get_size()):
                parts = MPI.WORLD.Recv(source=proc)
                for i in range(0,parts.shape[1]):
                    self._write_track(parts[:,i],s)
            self.add_times.append(time.time() - t0)
        else:
            save_ids = []
            for i in range(0,bunch.get_num_particles_local()):
                id = int(bunch.get_local_particles()[6,i])
                if (self.numer*id % self.denom) < self.numer:
                    save_ids.append(i)
            save_parts = numpy.zeros((7,len(save_ids)),'d')
            for i in range(0,len(save_ids)):
                save_parts[:,i] = bunch.get_local_particles()[:,save_ids[i]]
            MPI.WORLD.Send(save_parts,dest=0)
        self.count += 1
        if self.save_period:
            if (self.count % self.save_period == 0) and \
                (MPI.COMM_WORLD.Get_rank() == 0):
                t0 = time.time()
                self._copy()
                print "tracker auto copy time",time.time()-t0

    def _copy(self):
        if MPI.COMM_WORLD.Get_rank() == 0:
            t0 = time.time()
            destination = os.path.join(self.dest_dir,self.track_dir)
            if os.path.exists(destination):
                shutil.rmtree(destination)
            shutil.copytree(self.dir,destination)
            self.copy_time = time.time() - t0
            shutil.rmtree(self.base_dir)

    def close(self):
        self._copy()
        self.open = 0

    def show_statistics(self,filename=None):
        if MPI.COMM_WORLD.Get_rank() == 0:
            if filename == None:
                print "add times =",self.add_times
                print "copy time =",self.copy_time
            else:
                f = open(filename,'a')
                f.write('add times = %s\n' % string(self.add_times))
                f.write('copy times = %s\n' % string(self.copy_times))
                f.close()


def index2filename(index,dirname,suffix):
    i1 = index % 100
    i2 = ((index - i1) % 10000)/100
    i3 = ((index - i1 - i2*100) % 1000000)/10000
    i4 = ((index - i1 - i2*100 - i3*10000) % 100000000)/1000000
    return "%s/%02d/%02d/%02d/%02d.%s" % (dirname, i4, i3, i2, i1, suffix)
