#!/usr/bin/env python

from macro_bunch_store import Macro_bunch_store
import numpy
from mpi4py import MPI
from mpi4py import __version__ as mpi4py__version__
mpi4py_version = int(mpi4py__version__.split('.')[0])

import os.path
import tables
import time
import math

from synergia import loadfile_transpose, physics_constants, syn2_diagnostics 
import populate
import constraints

_need_random_init= True
_global_id_max = 0
class Macro_bunch:
    def __init__(self):  
        self.complete=0 
        self.particles = None
        self.local_num = None
        self.total_num = None # total number of macroparticles in the bunch
        self.total_current = None
        self.bunch_np=None  # total number of particles in the bunch
        self.units = None
        self.is_fixedz = None
        self.ref_particle = None
        self.periodic=None
        self.mass = None 
        self.charge = None 
        self.periodic=False
        self.periodic_z_size=None
        self.diagnostics=None
        self.bucket_num=None

# constructors     
    @classmethod
    def gaussian(klass, bunch_np,total_num,beam_parameters,diagnostics=None,bucket_num=0,periodic=False):
        bunch=klass() 
        bunch.init_gaussian(bunch_np,total_num,beam_parameters,periodic)        
        bunch.diagnostics=diagnostics
        bunch.bucket_num=bucket_num
        return bunch  
    
    @classmethod
    def gaussian_covariance(klass,bunch_np,total_num, beam_parameters,covariance,diagnostics=None,bucket_num=0,periodic=False):
        bunch=klass() 
        bunch.init_gaussian_covariance(bunch_np,total_num, beam_parameters,covariance,periodic)
        bunch.diagnostics= diagnostics 
        bunch.bucket_num=bucket_num     
        return bunch       
             
    @classmethod
    def from_bunch(klass, mbunch,diagnostics=None,bucket_num=0):
        bunch=klass() 
        bunch.init_from_bunch(mbunch)
        bunch.diagnostics=diagnostics  
        bunch.bucket_num=bucket_num     
        return bunch
        
    @classmethod
    def sphere(klass,num,radius,diagnostics=None,bucket_num=0):
        bunch=klass() 
        bunch.init_sphere(num,radius)
        bunch.diagnostics= diagnostics
        bunch.bucket_num=bucket_num        
        return bunch
            
    @classmethod
    def cylinder(klass,num,radius,length,diagnostics=None,bucket_num=0):
        bunch=klass() 
        bunch.init_cylinder(num,radius,length)
        bunch.diagnostics= diagnostics
        bunch.bucket_num=bucket_num         
        return bunch    
        
    @classmethod
    def test(klass,num_per_side,diagnostics=None,bucket_num=0,edge_length=1.0):
        bunch=klass() 
        bunch.init_test(num_per_side,edge_length)
        bunch.diagnostics= diagnostics
        bunch.bucket_num=bucket_num  
        return bunch
    
    
    @classmethod
    def test_am(klass,bunch_np,num_per_cell,shape,beam_parameters,diagnostics=None,bucket_num=0):
        bunch=klass() 
        bunch.init_test_am(bunch_np,num_per_cell,shape, beam_parameters)
        bunch.diagnostics= diagnostics
        bunch.bucket_num=bucket_num          
        return bunch
    
    # read in the bunch from a base file name and a bucket number.  The
    # actual file to read is constructed from the rank.
    @classmethod
    def read_bunch(self, filename, bucket_num=0, diagnostics=None):

        bunch = self()

        # I'll need my processor number
        myrank = MPI.COMM_WORLD.Get_rank()

        # get the file name
        hname = (filename + "-%03d-%05d.h5") % (bucket_num, myrank)

        # open the file
        h5file = tables.openFile(hname, "r")

        # read in the row that describes the bunch
        bdescr = [r for r in h5file.root.macro_bunch.iterrows()]

        # set the bunch characteristics
        bunch.bucket_num = bucket_num
        bunch.local_num = bdescr[0]['local_num']
        bunch.total_num = bdescr[0]['total_num']
        bunch.total_current = bdescr[0]['total_current']
        bunch.bunch_np = bdescr[0]['bunch_np']
        bunch.units = bdescr[0]['units']
        bunch.is_fixedz = bdescr[0]['is_fixedz']
        bunch.periodic = bdescr[0]['periodic']
        bunch.periodic_z_size = bdescr[0]['periodic_z_size']
        bunch.mass = bdescr[0]['mass']
        bunch.charge = bdescr[0]['charge']
        bunch.ref_particle = bdescr[0]['ref_particle']

        bunch.diagnostics = diagnostics

        # read in the particles
        bunch.particles = h5file.root.particles.read()

        # check that local_num is correct
        if bunch.particles.shape[1] != bunch.local_num:
            raise RuntimeError, "local_num %d  doesn't match the number of particles %d read on rank %d" % (bunch.local_num, bunch.particles.shape[1], myrank)

        # calculate and check total_num
        calc_total = MPI.COMM_WORLD.allreduce(bunch.local_num, op=MPI.SUM)
        if calc_total != bunch.total_num:
            raise RuntimeError, "total_num %d doesn't match computed total %d on rank %d" % (bunch.total_num, calc_total, myrank)

        h5file.close()

        return bunch

    def get_local_particles(self):
        return self.particles
    
    def get_num_particles_local(self):
        return self.local_num
    
    def get_bunch_np():
        return self.get_bunch_np()
    
    
    def get_store(self):
        if self.particles != None:
            return Macro_bunch_store(self.particles,
                                     self.local_num,
                                     self.total_num,
                                     self.mass,
                                     self.charge,
                                     self.total_current,
                                     self.bunch_np,
                                     self.units,
                                     self.ref_particle,
                                     self.is_fixedz)
        else:
            return None
        
        
    def get_longitudinal_period_size(self):   
        gamma = -self.get_store().get_ref_particle()[5]# beam frame
        # this gets the proper period for the fixedt frame when
        # the z coordinate of the particle store is in scaled units.
        return self.periodic_z_size*gamma*self.units[0]
        
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
        self.mass = 1.0
        self.charge = 1
        self.units = numpy.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = numpy.array([0.0,0.0,0.0,0.0,0.0,-1.1],'d')
        self.particles = numpy.zeros((7,local_num),'d')
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
	

    def init_test_am(self,bunch_np,num_per_cell,shape, beam_parameters):
        '''Particles uniformly distributed in a cube of size "size"'''
        cc=beam_parameters.get_covariances()
	
	 
        sigmax=math.sqrt(cc[0,0])
        sigmay=math.sqrt(cc[2,2])
        sigmaz=math.sqrt(cc[4,4])
        
	
        size = (2.*sigmax,2.*sigmay,2.*sigmaz)
        print" size=",size
        part_per_cell=int(num_per_cell)       
        (Cxy, Cxpyp, Cz, Czp) = beam_parameters.get_conversions()
        self.units = numpy.array([Cxy,Cxpyp,Cxy,Cxpyp,Cz,Czp],'d')
        offset = (beam_parameters.offset_x_m*Cxy,beam_parameters.offset_y_m*Cxy,beam_parameters.offset_z*Cz)
        num = (part_per_cell*shape[0],part_per_cell*shape[1],part_per_cell*shape[2])
        local_num = num[0]*num[1]*num[2]
        total_num = local_num # fix me for mpi
        self.mass = beam_parameters.get_mass()
        self.charge = beam_parameters.get_charge()
        self.is_fixedz = 1
        self.total_current = bunch_np*physics_constants.PH_MKS_e *beam_parameters.scaling_frequency_Hz
        self.bunch_np=bunch_np
        self.ref_particle = numpy.zeros((6,),'d')
        self.ref_particle[5] = -beam_parameters.get_gamma()        
        self.particles = numpy.zeros((7,local_num),'d')
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
  

    def init_sphere(self,num,radius):
        '''Particles uniformly distributed in a sphere of radius "radius"'''
        numpy.random.RandomState([17+MPI.COMM_WORLD.Get_rank()*51,59+MPI.COMM_WORLD.Get_rank()*23])
        offset = (0.0,0.0,0.0)
        local_num = num/MPI.COMM_WORLD.Get_size() #jfa: could be more precise...
        total_num = MPI.COMM_WORLD.Allreduce(local_num,MPI.SUM)
        total_current = 1.0
        self.mass=1.0
        self.charge=1
        self.units = numpy.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = numpy.array([0.0,0.0,0.0,0.0,0.0,-1.1],'d')
        self.particles = numpy.zeros((7,local_num),'d')
        self.is_fixedz = 1
        index = 0
        added = 0
        discarded = 0
        chunk_size = 1000
        t0 = time.time()
        p = (numpy.random.random([6,chunk_size])-0.5)*2.0*radius
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
                p = (numpy.random.random([6,chunk_size])-0.5)*2.0*radius
                index = 0
        t1 = time.time()
#        print "pi =",6.0*added/(1.0*added+discarded),"in",t1-t0,"secs"
        self.local_num = local_num
        self.total_num = total_num
        self.total_current = 1.0
        self.complete = 1

    def init_cylinder(self,num,radius,length):
        '''Particles uniformly distributed in a cylinder of radius "radius"
        and length 2pi'''
        offset = (0.0,0.0,0.0)
        local_num = num
        total_num = local_num # fix me for mpi
        total_current = 1.0
        self.mass=1.0
        self.charge=1
        self.units = numpy.array([1.0,1.0,1.0,1.0,1.0,1.0],'d')
        self.ref_particle = numpy.array([0.0,0.0,0.0,0.0,0.0,-1.1],'d')
        self.particles = numpy.zeros((7,local_num),'d')
        self.is_fixedz = 1
        index = 0
        added = 0
        discarded = 0
        chunk_size = 1000
        t0 = time.time()
        p = (numpy.random.random([6,chunk_size])-0.5)*2.0*radius
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
                p = (numpy.random.random([6,chunk_size])-0.5)*2.0*radius
                index = 0
        t1 = time.time()
        self.local_num = local_num
        self.total_num = total_num
        self.total_current = total_current
        self.complete = 1

    def _split_num_particles(self, total_num, divisions):
        div_num = total_num/divisions
        remainder = total_num % divisions
        nums = []
        offsets = []
        offset = 0 
        for i in range(0,divisions):
            if i < remainder:
                num = div_num + 1
            else:
                num = div_num
            nums.append(num)
            offsets.append(offset)
            offset += num
        return nums,offsets

    
    def init_gaussian(self,bunch_np,total_num,beam_parameters, periodic=False):
        
        if periodic:
            if (beam_parameters.get_z_length()==None):
                raise RuntimeError, " z_length is needed for a periodic bunch!"
            self.periodic_z_size=beam_parameters.get_z_length()
        self.periodic=periodic   
            
        (Cxy, Cxpyp, Cz, Czp) = beam_parameters.get_conversions()
        self.units = numpy.array([Cxy,Cxpyp,Cxy,Cxpyp,Cz,Czp],'d')
        self.total_num = total_num       
        nums, offsets = self._split_num_particles(total_num,MPI.COMM_WORLD.Get_size())
        self.local_num = nums[MPI.COMM_WORLD.Get_rank()]
        global _global_id_max
        id_offset = _global_id_max + offsets[MPI.COMM_WORLD.Get_rank()]
        _global_id_max += total_num
        self.mass = beam_parameters.get_mass()
        self.charge = beam_parameters.get_charge()
        self.bunch_np=bunch_np

        self.is_fixedz = 1
        self.ref_particle = numpy.zeros((6,),'d')
        self.ref_particle[5] = -beam_parameters.get_gamma()
        self.particles = numpy.zeros((7,self.local_num),'d')
        # jfa: the following parallel seed algorithm is ad hoc
        #      Need to get better algorithm, probably from SPRNG
        seed_offset = int(time.time())
        long_seed = (1000+5*(MPI.COMM_WORLD.Get_rank()+seed_offset))*((MPI.COMM_WORLD.Get_rank()+seed_offset)+7)-1
        seed = int(long_seed % 2**32)
        global _need_random_init
        if beam_parameters.get_transverse():
            t_length=beam_parameters.get_t_length()    
            populate.populate_transverse_gaussian(self.particles,
                beam_parameters.get_means(),
                beam_parameters.get_covariances(),t_length,id_offset,seed,
                _need_random_init)
            _need_random_init= False 
            # current below should be used for a transverse beam, and is needed ONLY in the gaussian solver
            self.total_current = bunch_np*physics_constants.PH_MKS_e*beam_parameters.get_beta()* \
                  physics_constants.PH_MKS_c/beam_parameters.get_z_length()       
        else:
            if beam_parameters.get_z_peaks() == 1:
                populate.populate_6d_gaussian(self.particles,
                    beam_parameters.get_means(),
                    beam_parameters.get_covariances(),id_offset,seed,
                    _need_random_init)
                _need_random_init= False
                # you shouldn't  use anywhere the total_current defined below
                self.total_current = bunch_np*physics_constants.PH_MKS_e *beam_parameters.scaling_frequency_Hz
               # self.completed=1 
            else:
               # raise RuntimeError, "please modify the init of a bunch with multiple z_peaks"
                #delta_z = 2.0*math.pi
                if (beam_parameters.get_z_length()==None):
                    raise RuntimeError, " z_length is needed for multiple peaks!" 
                delta_z =beam_parameters.get_t_length()
                num_peaks = beam_parameters.get_z_peaks()
                old_means = beam_parameters.get_means()
                new_means = numpy.zeros([len(old_means)],'d')
                new_means[:] = old_means[:]
                peak_nums, peak_offsets = \
                    self._split_num_particles(self.local_num,num_peaks)
                peak_offsets.append(self.local_num)
                for peak in range(0,num_peaks):
                    new_means[4] = old_means[4] + \
                        ((2*peak+1)/(2.0*num_peaks) - 0.5)*delta_z
                    tmp = numpy.zeros((7,peak_nums[peak]),'d')
                    populate.populate_6d_gaussian(tmp,new_means,
                        beam_parameters.get_covariances(),
                        id_offset+peak_offsets[peak],seed,
                        _need_random_init)
                    self.particles[:,peak_offsets[peak]:peak_offsets[peak+1]] = tmp
                    _need_random_init= False
                self.total_current = bunch_np*physics_constants.PH_MKS_e*beam_parameters.get_beta()* \
                  physics_constants.PH_MKS_c/beam_parameters.get_z_length()   
            if periodic:
                t_length=beam_parameters.get_t_length() 
                constraints.apply_longitudinal_periodicity(self.get_store(),t_length)
                    
    def init_gaussian_covariance(self,bunch_np,total_num, beam_parameters,covariance,periodic=False):
        
        if periodic:
            if (beam_parameters.get_z_length()==None):
                raise RuntimeError, " z_length is needed for a periodic bunch!"
            self.periodic_z_size=beam_parameters.get_z_length()
        self.periodic=periodic
            
        (Cxy, Cxpyp, Cz, Czp) = beam_parameters.get_conversions()
        self.units = numpy.array([Cxy,Cxpyp,Cxy,Cxpyp,Cz,Czp],'d')
        self.total_num = total_num
        nums, offsets = self._split_num_particles(total_num,MPI.COMM_WORLD.Get_size())
        self.local_num = nums[MPI.COMM_WORLD.Get_rank()]
        global _global_id_max
        id_offset = _global_id_max + offsets[MPI.COMM_WORLD.Get_rank()]
        _global_id_max += total_num
        self.mass = beam_parameters.get_mass()
        self.charge = beam_parameters.get_charge()
        self.bunch_np=bunch_np
        self.total_current = bunch_np*physics_constants.PH_MKS_e *beam_parameters.scaling_frequency_Hz
        self.periodic=periodic
        self.is_fixedz = 1
        self.ref_particle = numpy.zeros((6,),'d')
        self.ref_particle[5] = -beam_parameters.get_gamma()
        self.particles = numpy.zeros((7,self.local_num),'d')
        seed_offset = int(time.time())
        long_seed = (1000+5*(MPI.COMM_WORLD.Get_rank()+seed_offset))*((MPI.COMM_WORLD.Get_rank()+seed_offset)+7)-1
        seed = int(long_seed % 2**32)
        global _need_random_init
        if beam_parameters.get_transverse():
            raise RuntimeError, "init_gaussian_covariance does not work with transverse beams"
        else:
            if beam_parameters.get_z_peaks() == 1:
                populate.populate_6d_gaussian(self.particles,
                    beam_parameters.get_means(),
                    covariance,id_offset,seed,
                    _need_random_init)
                _need_random_init= False
            else:
                raise RuntimeError, "init_gaussian_covariance does not work with z_peaks <> 1"
            if periodic:
                t_length=beam_parameters.get_t_length() 
                constraints.apply_longitudinal_periodicity(self.get_store(),t_length)
            
    def inject(self, bunch):
        if (self.is_fixedz != bunch.is_fixedz):
            raise RuntimeError, "injected bunch must have same fixedz status as parent"
        if (self.ref_particle != bunch.ref_particle):
            raise RuntimeError, "injected bunch must have same reference particle as parent"
        
        new_total_num = self.total_num + bunch.total_num
        new_local_num = self.local_num + bunch.local_num
        new_total_current = self.total_current + bunch.total_current
        new_particles = numpy.zeros((7,new_local_num),'d')
        new_particles[:,0:self.local_num] = self.particles
        new_particles[:,self.local_num:new_local_num] = bunch.particles
        self.particles = new_particles
        self.total_num = new_total_num
        self.local_num = new_local_num
        self.total_current = new_total_current
        
    def init_from_bunch(self, bunch):
        self.units = bunch.units
        self.mass = bunch.mass
        self.charge = bunch.charge
        self.particles = bunch.get_local_particles()
        self.local_num = bunch.get_num_particles_local()
        self.total_num = bunch.total_num
        self.total_current = bunch.total_current
        self.bunch_np=bunch.get_bunch_np()
        self.ref_particle = bunch.ref_particle
        self.is_fixedz = 1
        
    def write_particles(self,filename,compress_level=1):
        old_pytables = False
        try:
            if tables.__version__.split('.')[0] == '1':
                old_pytables = True
        except:
            pass        
        if MPI.COMM_WORLD.Get_rank() == 0:
            h5filename = os.path.splitext(filename)[0] + '.h5'
            f = tables.openFile(h5filename,mode = "w")
            filter = tables.Filters(complevel=compress_level)
            if old_pytables:
                atom = tables.Atom(dtype='Float64',shape=(7,0))
                earray = f.createEArray(f.root,'particles',atom,'Float64',
                                        filters = filter)
            else:
                atom = tables.Float64Atom()
                earray = f.createEArray(f.root,'particles',atom,(7,0),
                                        filters = filter)
            particles = self.particles
            earray.append(particles)
            for proc in xrange(1,MPI.COMM_WORLD.Get_size()):
                if mpi4py_version > 0:
                    parts = MPI.COMM_WORLD.recv(source=proc)
                else:
                    parts = MPI.COMM_WORLD.Recv(source=proc)
                if parts.shape[1] > 0:
                    earray.append(parts)
            f.close()
        else:
            if mpi4py_version > 0:
                MPI.COMM_WORLD.send(self.particles,dest=0)
            else:
                MPI.COMM_WORLD.Send(self.particles,dest=0)
          

    def write_particles_text(self,filename):
        if MPI.COMM_WORLD.Get_rank() == 0:
            f = open(filename,"w")
            for proc in xrange(1,MPI.COMM_WORLD.Get_size()):
                if mpi4py_version > 0:
                    parts = MPI.COMM_WORLD.recv(source=proc)
                else:
                    parts = MPI.COMM_WORLD.Recv(source=proc)
                for i in range(0,parts.shape[1]):
                    f.write("%g %g %g %g %g %g %g\n" % \
                            tuple(parts[:,i]))
            parts = self.particles
            for i in range(0,parts.shape[1]):
                f.write("%g %g %g %g %g %g %g\n" % \
                        tuple(parts[:,i]))
            f.close()
        else:
            if mpi4py_version > 0:
                MPI.COMM_WORLD.send(self.particles,dest=0)
            else:
                MPI.COMM_WORLD.Send(self.particles,dest=0)

    def read_particles(self,filename,total_current, beam_parameters):
        (Cxy, Cxpyp, Cz, Czp) = beam_parameters.get_conversions()
        self.units = numpy.array([Cxy,Cxpyp,Cxy,Cxpyp,Cz,Czp],'d')
        self.total_current = total_current
        self.is_fixedz = 1
        self.ref_particle = numpy.zeros((6,),'d')
        self.ref_particle[5] = -beam_parameters.get_gamma()
        self.particles = loadfile_transpose(filename)
        self.local_num = self.particles.shape[1]
        self.total_num = self.local_num
        print "read",self.local_num,"particles"
    
    def add_diagnostics(self,s):
        self.diagnostics.add(s,self)

    
    # write the current processors notion of a bunch including
    # particles and characteristics to an hdf5 file to be read in
    # later.  filename is used as the template for the actual
    # name which will be filename-bbb-nnnnn.h5 bbb is the bucket
    # number, nnnnnn is the rank.
    def write_bunch(self, filename):
        # I'll need my processor number
        myrank = MPI.COMM_WORLD.Get_rank()

        # create the hdf5 file
        hname = (filename + "-%03d-%05d.h5") % (self.bucket_num, myrank)
        h5file = tables.openFile(hname, "w",
                                 title="rank %d bunch bucket %d" %
                                        (myrank, self.bucket_num))

        # this class mirrors the Macro bunch definition
        class bunch_prop(tables.IsDescription):
            local_num = tables.Int64Col()
            total_num = tables.Int64Col()
            total_current = tables.Float64Col()
            bunch_np = tables.Float64Col()
            bucket_num = tables.Int64Col()
            units = tables.Float64Col(6)
            is_fixedz = tables.BoolCol()
            ref_particle = tables.Float64Col(6)
            periodic = tables.BoolCol()
            periodic_z_size = tables.BoolCol()
            mass = tables.Float64Col()
            charge = tables.Int64Col()

        compress_level=1
        filter = tables.Filters(complevel=compress_level)

        # make a table (overkill) for the bunch characteristics
        btab = h5file.createTable(h5file.root, 'macro_bunch', bunch_prop, 'bunch characteristics', filters=filter)

        b=btab.row
        b['local_num'] = self.local_num
        b['total_num'] = self.total_num
        b['total_current'] = self.total_current
        b['bunch_np'] = self.bunch_np
        b['bucket_num'] = self.bucket_num
        b['is_fixedz'] = self.is_fixedz
        b['units'] = self.units
        b['ref_particle'] = self.ref_particle
        b['periodic'] = self.periodic
        b['periodic_z_size'] = self.periodic_z_size
        b['mass'] = self.mass
        b['charge'] = self.charge
        b.append()

        btab.flush()

        # write out the particles
        p = h5file.createArray(h5file.root, 'particles', self.particles)
        h5file.close()

class Multiple_bunches:
    def __init__(self, bunches,  bunch_spacing):
        self.bunches= bunches     
        self.bunch_spacing=bunch_spacing
        self.num_buckets=len(bunches)
        
      
