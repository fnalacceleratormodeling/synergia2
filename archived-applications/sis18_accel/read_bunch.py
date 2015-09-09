import os
import synergia
import numpy as np
import tables
from mpi4py import MPI

# load the particles that will be used for the simulation
# The particles file is a text file with particle coordinates
# defined with the MAD-X conventions: X PX Y PY T PT
# Read this in using numpy's loadtxt command
# particle coordinates are converted to Synergia conventions

# input arguments:
#    particles_file: the file name
#    reference particle: the lattice reference particle for kinematic conversions
#    real_particles: the real charge of the bunch
#    bucket_length: the longitudinal length of the bucket
#    comm: the Commxx communicator object for this bunch
#    verbose: be chatty about what's happening
#  
def read_bunch(particles_file, refpart, real_particles, bucket_length, comm, verbose=False):
    name,extension = os.path.splitext(particles_file)
    if extension == ".mxtxt":
        return read_txt_particles(particles_file, refpart, real_particles, bucket_length, comm, True, verbose)
    elif extension == ".txt":
        return read_txt_particles(particles_file, refpart, real_particles, bucket_length, comm, False, verbose)
    elif extension == ".h5":
        return read_h5_particles(particles_file, refpart, real_particles, bucket_length, comm, verbose)
    else:
        raise RuntimeError, "unrecognized file format: %s"%extension

#====================================================================

# if madx_format is True, the particles are in madX units, otherwise they are in
# synergia units

def read_txt_particles(particles_file, refpart, real_particles, bucket_length, comm, madx_format, verbose):
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0:
        if madx_format:
            print "Loading madX particles from txt file: ", particles_file
        else:
            print "Loading Synergia particles from txt file: ", particles_file

    if myrank == 0:
        particles = np.loadtxt(particles_file)
        num_total_particles = particles.shape[0]
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particles.shape[1] != 6) and (particles.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)
        
        
        if madx_format:
            # numpy manipulations to convert kinematics
            # convert MAD-X T=-c*dt to Synergia c*ct
            particles[:,4] = -particles[:,4]
            # convert MAD-X Delta-E/pc to Synergia delta-p/p
            # sqrt(((dE/p0c)+(E0/p0c))**2 - (m/p0c)**2) - (p0c/p0c)
            m_over_pc = pmass/p0c
            E_0_over_pc = E_0/p0c
            particles[:,5] = np.sqrt( (particles[:,5] + E_0_over_pc) *
                                      (particles[:,5] + E_0_over_pc) - m_over_pc**2 ) - 1.0
        

        # if there are no IDs, append particle ID column
        if particles.shape[1] != 7:
            particles_w_id = np.column_stack((particles,
                                              np.arange(num_total_particles, dtype='d')))
        else:
            particles_w_id = particles
            
            if myrank == 0:
                print "Read ", num_total_particles, " particles"


    # create a bunch with the correct number of macro particles
    bunch = synergia.bunch.Bunch(
        refpart,
        num_total_particles, real_particles, comm,
        bucket_length)

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particles_w_id[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles_w_id[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]
    return bunch

#==========================================================

def read_h5_particles(particles_file, refpart, real_particles, bucket_length, comm, verbose):
    four_momentum = refpart.get_four_momentum()
    pmass = four_momentum.get_mass()
    E_0 = four_momentum.get_total_energy()
    p0c = four_momentum.get_momentum()

    myrank = comm.get_rank()
    mpisize = comm.get_size()
    
    if myrank==0:
        print "Loading particles from h5 file: ", particles_file

    if myrank == 0:
        h5 = tables.openFile(particles_file)
        # use explicit int conversion otherwise there seems to
        # be a typepython->C++ type  mismatch of numpy.int64->int
        num_total_particles = int(h5.root.particles.shape[0])
        if verbose:
            print "Total of  ", num_total_particles, " particles from file"
        # broadcast num particles to all nodes
        MPI.COMM_WORLD.bcast(num_total_particles, root=0)
    else:
        num_total_particles = None
        num_total_particles = MPI.COMM_WORLD.bcast(num_total_particles, root=0)

    if myrank == 0:
        particles = h5.root.particles
        # make sure the data has the correct shape, either [n,6] without
        # particles IDs or [n,7] with particle IDs.
        if (particles.shape[1] != 7):
            raise RuntimeError, "input data shape %shas incorrect number of particle coordinates"%repr(particles.shape)


    # create a bunch with the correct number of macro particles
    bunch = synergia.bunch.Bunch(
        refpart,
        num_total_particles, real_particles, comm,
        bucket_length, 0)

    local_num = bunch.get_local_num()
    local_particles = bunch.get_local_particles()

    # Each processor will have a possibly different number of local particles.
    # rank 0 has to find out how many each of them has and distribute them
    n_particles_by_proc = MPI.COMM_WORLD.gather(local_num, 0)
    if myrank == 0:
        # copy in my particles
        this_rank_start = 0
        local_particles[:,:] = particles[0:local_num, :]
        this_rank_start += local_num
        # send particles out to other ranks
        for r in range(1, mpisize):
            this_rank_end = this_rank_start+n_particles_by_proc[r]
            MPI.COMM_WORLD.send(obj=particles[this_rank_start:this_rank_end, :],
                                dest=r)
            this_rank_start += n_particles_by_proc[r]
    else:
        # I'm not rank 0.  Receive my particles
        lp = MPI.COMM_WORLD.recv(source=0)
        local_particles[:,:] = lp[:,:]

    return bunch

#================================================================

def print_bunch_stats(bunch):
    coord_names = ("x", "xp", "y", "yp", "c*dt", "dp/p")

    myrank = bunch.get_comm().get_rank()
    means = synergia.bunch.Core_diagnostics().calculate_mean(bunch)
    stds = synergia.bunch.Core_diagnostics().calculate_std(bunch, means)
    if myrank == 0:
        print "%20s   %20s   %20s"%("coord","mean","rms")
        print "%20s   %20s   %20s"%("====================",
                                    "====================",
                                    "====================")
        for i in range(6):
            print "%20s   %20.12e   %20.12e"%(coord_names[i], means[i], stds[i])

#=========================================================
