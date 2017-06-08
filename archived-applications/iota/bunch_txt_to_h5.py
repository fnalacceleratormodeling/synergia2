#!/usr/bin/env python
# -*- coding: utf-8 -*-
import os
import sys
import synergia
import numpy
from mpi4py import MPI
#from modes_options import opts


#def read_txt_particles(particles_file):
   #particles = numpy.loadtxt(particles_file)
   #if (particles.shape[1] != 7):
        #raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)

   #num_total_particles = particles.shape[0]
   #print "num_total_particles=",num_total_particles
   ## fake ref particle
   #refpart =synergia.foundation.Reference_particle(1, 1., 5.);
   #real_particles=100.;
   #comm = synergia.utils.Commxx(True)
   #bunch = synergia.bunch.Bunch(refpart,  num_total_particles, real_particles, comm)
   
   
  # Commxx_sptr commx=comms[i]; 
  #       Bunch_sptr bunch_sptr=Bunch_sptr(new Bunch(stepper_sptr->get_lattice_simulator().get_lattice().get_reference_particle(),
  #           opts.num_macroparticles, opts.num_real_particles, commx));  

#def read_txt_particles(particles_file, refpart, real_particles, bucket_length, comm, madx_format, verbose):
    #four_momentum = refpart.get_four_momentum()
    #pmass = four_momentum.get_mass()
    #E_0 = four_momentum.get_total_energy()
    #p0c = four_momentum.get_momentum()

    #myrank = comm.get_rank()
    #mpisize = comm.get_size()
    
    #if myrank==0:
        #if madx_format:
            #print "Loading madX particles from txt file: ", particles_file
        #else:
            #print "Loading Synergia particles from txt file: ", particles_file

    #if myrank == 0:
        #particles = np.loadtxt(particles_file)
        #num_total_particles = particles.shape[0]
        ## broadcast num particles to all nodes
        #MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    #else:
        #num_total_particles = None
        #num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    #if myrank == 0:
        ## make sure the data has the correct shape, either [n,6] without
        ## particles IDs or [n,7] with particle IDs.
        #if (particles.shape[1] != 6) and (particles.shape[1] != 7):
            #raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)
        
        
        #if madx_format:
            ## numpy manipulations to convert kinematics
            ## convert MAD-X T=-c*dt to Synergia c*ct
            #particles[:,4] = -particles[:,4]
            ## convert MAD-X Delta-E/pc to Synergia delta-p/p
            ## sqrt(((dE/p0c)+(E0/p0c))**2 - (m/p0c)**2) - (p0c/p0c)
            #m_over_pc = pmass/p0c
            #E_0_over_pc = E_0/p0c
            #particles[:,5] = np.sqrt( (particles[:,5] + E_0_over_pc) *
                                      #(particles[:,5] + E_0_over_pc) - m_over_pc**2 ) - 1.0
        

        ## if there are no IDs, append particle ID column
        #if particles.shape[1] != 7:
            #particles_w_id = np.column_stack((particles,
                                              #np.arange(num_total_particles, dtype='d')))
        #else:
            #particles_w_id = particles
            
            #if myrank == 0:
                #print "Read ", num_total_particles, " particles"


    ## create a bunch with the correct number of macro particles
    ##bunch = synergia.bunch.Bunch(
    ##    refpart,
    ##    num_total_particles, real_particles, comm)
    ##bunch.set_z_period_length(bucket_length)
    
    ##Note: Synergia bunch constructor updated - commit 077b99d7 - 11/17/2016
    ##Using old constructor throws an ArgumentError of a non-standard type.
    ## Using a try and except to handle both instances.
    #try:
        ## try the original constructor
        #bunch = synergia.bunch.Bunch(
            #refpart,
            #num_total_particles, real_particles, comm,
            #bucket_length)
    #except Exception, e:
        ##look to see if it's an ArgumentError by evaluating the traceback
        #if (not str(e).startswith("Python argument types in")):
            #raise
        #else:
            ## use the new constructor
            #print "Using updated bunch constructor"
            #bunch = synergia.bunch.Bunch(
                #refpart,
                #num_total_particles, real_particles, comm)
            #bunch.set_z_period_length(bucket_length)

    #local_num = bunch.get_local_num()
    #local_particles = bunch.get_local_particles()

    ## Each processor will have a possibly different number of local particles.
    ## rank 0 has to find out how many each of them has and distribute them
    #n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    #if myrank == 0:
        ## copy in my particles
        #this_rank_start = 0
        #local_particles[:,:] = particles_w_id[0:local_num, :]
        #this_rank_start += local_num
        ## send particles out to other ranks
        #for r in range(1, mpisize):
            #this_rank_end = this_rank_start+n_particles_by_proc[r]
            #MPI.COMM_WORLD.send(obj=particles_w_id[this_rank_start:this_rank_end, :],
                                #dest=r)
            #this_rank_start += n_particles_by_proc[r]
    #else:
        ## I'm not rank 0.  Receive my particles
        #lp = MPI.COMM_WORLD.recv(source=0)
        #local_particles[:,:] = lp[:,:]
    #return bunch

#==========================================================

if __name__ == "__main__":

    if len(sys.argv) < 2:
       print " nput the txt file name after the comand line "
       sys.exit(1)

    txtfile = sys.argv[1]

    print "bunch txt file is: ",txtfile
    particles = numpy.loadtxt(txtfile)
    if (particles.shape[1] != 7):
        raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)

    num_total_particles = particles.shape[0]
    print "num_total_particles=",num_total_particles
    # fake ref particle, it should not matter
    refpart =synergia.foundation.Reference_particle(1, 1., 5.);
    real_particles=100.;
    comm = synergia.utils.Commxx(True)
    bunch = synergia.bunch.Bunch(refpart,  num_total_particles, real_particles, comm)
    local_particles = bunch.get_local_particles()
    print "local_particles shape=",local_particles.shape[0], local_particles.shape[1]
    print "particles shape=",particles.shape[0], particles.shape[1]
    local_particles[:,:]=particles[:,:]
    
    outputfile=txtfile[0:-4]+".h5"
    print "outputfile=",outputfile
    particles_diagnostics=synergia.bunch.Diagnostics_particles(outputfile)
    particles_diagnostics.set_bunch(bunch)
    particles_diagnostics.update_and_write()
            
    
  