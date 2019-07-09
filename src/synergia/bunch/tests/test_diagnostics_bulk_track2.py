import sys
import synergia

from mpi4py import MPI
from synergia.foundation import Reference_particle, Four_momentum, pconstants
from synergia.bunch import Bunch
from synergia.bunch import Diagnostics_bulk_track
from synergia.utils import Commxx, Hdf5_file, decompose_1d
import synergia.convertors
from synergia.lattice import Lattice, Lattice_element
from synergia.simulation import Independent_stepper, Bunch_simulator, Propagator
import numpy
from nose.tools import *

mass = pconstants.mp
momentum = 1.0
total_num = 100
real_num = 1.0e9
proton_charge = 1


comm = Commxx()

def test_tracks1():
    four_momentum = Four_momentum(mass)
    four_momentum.set_momentum(momentum)
    reference_particle = Reference_particle(proton_charge, four_momentum)
    bunch = Bunch(reference_particle, total_num, real_num, comm)
    local_particles = bunch.get_local_particles()
    local_num = bunch.get_local_num()

    lattice = Lattice("foo")
    d = Lattice_element("drift", "d")
    d.set_double_attribute("l", 1.0)
    lattice.append(d)
    lattice.set_reference_particle(reference_particle)

    stepper = Independent_stepper(lattice, 1, 1)

    bunch_simulator = Bunch_simulator(bunch)

    bunch_simulator.add_per_turn(Diagnostics_bulk_track("dummy.h5", total_num))

    # set particles to flare out
    # particle 0 goes straight
    # particle 1 goes out 1 mm/turn
    # particle 2 goes out 2 mm/turn
    #  ...

    l = lattice.get_length()
    dx = 0.001
    (offsets, counts) = decompose_1d(comm, total_num)
    mystart = offsets[comm.get_rank()]
    mycount = counts[comm.get_rank()]
    assert_equals(mycount, local_num, msg="mycount != local_num")
    for pnum in range(local_num):
        i = mystart+pnum
        local_particles[pnum, 0:6] = 0.0
        local_particles[pnum, 1] = (i*dx) / numpy.sqrt(l**2 + (i*dx)**2)

    propagator = Propagator(stepper)

    propagator.propagate(bunch_simulator, 100, 0, 0)

    if (comm.get_rank() != 0): return

    h5 = Hdf5_file("dummy.h5", Hdf5_file.read_only)
    dims = h5.get_dims("track_coords")
    assert_equals(dims[0], 100, "number of particles")
    assert_equals(dims[1], 7, "dimensions")
    assert_equals(dims[2], 101, "number turns + 1")

    npart = dims[0]
    nturns = dims[2]
    track_coords = h5.read_array3d("track_coords")

    for p in range(npart):
        for t in range(nturns):
            # x coordinate for particle p increases 0.001*p per turn
            assert_almost_equal(track_coords[p, 0, t], 0.001*p*t, places=10, msg="x coordinate particle %d turn %d"%(p, t))

    # delete stuff to close diagnostics files
    del propagator
    del stepper
    del bunch_simulator
    del lattice


#------------------------------------------------------------------------------------------------    

def test_tracks_lost_particles():
    four_momentum = Four_momentum(mass)
    four_momentum.set_momentum(momentum)
    reference_particle = Reference_particle(proton_charge, four_momentum)
    bunch = Bunch(reference_particle, total_num, real_num, comm)
    local_particles = bunch.get_local_particles()
    local_num = bunch.get_local_num()

    (offsets, counts) = decompose_1d(comm, total_num)
    mystart = offsets[comm.get_rank()]
    mycount = counts[comm.get_rank()]
    assert_equals(mycount, local_num, msg="mycount != local_num")

    lattice = Lattice("foo")
    d = Lattice_element("drift", "d")
    d.set_double_attribute("l", 1.0)
    lattice.append(d)
    mrk = Lattice_element("marker", "mask")
    mrk.set_string_attribute("aperture_type", "circular")
    mrk.set_double_attribute("circular_aperture_radius", 0.1001)
    lattice.append(mrk)

    lattice.set_reference_particle(reference_particle)

    stepper = Independent_stepper(lattice, 1, 1)

    bunch_simulator = Bunch_simulator(bunch)

    bunch_simulator.add_per_turn(Diagnostics_bulk_track("dummy.h5", total_num))

    # set particles to increase 1 mm/turn
    # start particles at increasing offsets
    # aperture at end of drift will start to cut particles

    l = lattice.get_length()
    dx = 0.001
    npx = dx/numpy.sqrt(l**2 + dx**2)

    for pnum in range(local_num):
        i = mystart+pnum
        local_particles[pnum, 0:6] = 0.0
        local_particles[pnum, 0] = (total_num-i-1)*0.001
        local_particles[pnum, 1] = npx

    propagator = Propagator(stepper)

    propagator.propagate(bunch_simulator, 100, 0, 0)

    if (comm.get_rank() == 0): return

    h5 = Hdf5_file("dummy.h5", Hdf5_file.read_only)
    dims = h5.get_dims("track_coords")
    assert_equals(dims[0], 100, "number of particles")
    assert_equals(dims[1], 7, "dimensions")
    assert_equals(dims[2], 101, "number turns + 1")

    npart = dims[0]
    nturns = dims[2]
    track_coords = h5.read_array3d("track_coords")

    for p in range(npart):
        # x coordinate for particle p starts at (99-p)*0.001 and increases 1 mm/turn until it goes over 0.1001
        for t in range(nturns):
            x = (99-p)*0.001 + 0.001*t
            if x < 0.1001:
                assert_almost_equal(track_coords[p, 0, t], x, places=10, msg="x coordinate particle %d turn %d"%(p, t))
            else:
                assert_equals(track_coords[p, 0, t], 0.0, "x coordinate particle %d after being cut"%p)


    # delete stuff to close diagnostics files
    del propagator
    del stepper
    del bunch_simulator
    del lattice

#----------------------------------------------------------------------------------------------

def test_tracks_lost_particles_with_offset():
    four_momentum = Four_momentum(mass)
    four_momentum.set_momentum(momentum)
    reference_particle = Reference_particle(proton_charge, four_momentum)
    bunch = Bunch(reference_particle, total_num, real_num, comm)
    local_particles = bunch.get_local_particles()
    local_num = bunch.get_local_num()

    lattice = Lattice("foo")
    d = Lattice_element("drift", "d")
    d.set_double_attribute("l", 1.0)
    lattice.append(d)
    mrk = Lattice_element("marker", "mask")
    mrk.set_string_attribute("aperture_type", "circular")
    mrk.set_double_attribute("circular_aperture_radius", 0.1001)
    lattice.append(mrk)

    lattice.set_reference_particle(reference_particle)

    stepper = Independent_stepper(lattice, 1, 1)

    bunch_simulator = Bunch_simulator(bunch)

    # save the last 20 tracks
    bunch_simulator.add_per_turn(Diagnostics_bulk_track("dummy.h5", 20, total_num-20))

    # set particles to increase 1 mm/turn
    # start particles at increasing offsets
    # aperture at end of drift will start to cut particles

    l = lattice.get_length()
    dx = 0.001
    npx = dx/numpy.sqrt(l**2 + dx**2)

    myrank = comm.get_rank()
    size = comm.get_size()
    (local_starts, local_nums) = decompose_1d(comm, total_num)
    (offset_offsets, offset_counts) = decompose_1d(comm, total_num-20)
    (xxxx, try_counts) = decompose_1d(comm, 20)

    starts = [local_starts[k]+offset_counts[k] for k in range(size)]
    counts = list(try_counts)
    for k in range(size):
        if starts[k]+try_counts[k] > local_starts[k]+local_nums[k]:
            counts[k] = local_nums[k] - offset_counts[k]
            starts[k] = local_starts[k] + offset_counts[k]

    mystart = starts[myrank]
    mycount = counts[myrank]

    # make list of particle IDs selected with offset and track count
    pids = []
    for k in range(size):
        pids.extend(list(range(starts[k], starts[k]+counts[k])))

    # now filling local particle array for all particles
    for pnum in range(local_num):
        i = local_starts[myrank]+pnum
        local_particles[pnum, 0:6] = 0.0
        local_particles[pnum, 0] = (total_num-i-1)*0.001
        local_particles[pnum, 1] = npx

    propagator = Propagator(stepper)

    propagator.propagate(bunch_simulator, 100, 0, 0)

    if (comm.get_rank() == 0): return

    h5 = Hdf5_file("dummy.h5", Hdf5_file.read_only)
    dims = h5.get_dims("track_coords")
    assert_equals(dims[0], 20, "number of particles")
    assert_equals(dims[1], 7, "dimensions")
    assert_equals(dims[2], 101, "number turns + 1")

    npart = dims[0]
    nturns = dims[2]
    track_coords = h5.read_array3d("track_coords")

    # because of the way particles get distributed, there might not be
    # npart real particles in the array

    for pnum in range(len(pids)):
        # saving the last 20 particles of the previous test so
        # x coordinate for particle p starts at (99-p)*0.001 and increases 1 mm/turn until it goes over 0.1001
        for t in range(nturns):
            p = pids[pnum]
            x = (99-p)*0.001 + 0.001*t
            if x < 0.1001:
                assert_almost_equal(track_coords[pnum, 0, t], x, places=10, msg="x coordinate particle %d turn %d"%(p, t))
            else:
                assert_equals(track_coords[pnum, 0, t], 0.0, "x coordinate particle %d after being cut"%p)

    # delete stuff to close diagnostics files
    del propagator
    del stepper
    del bunch_simulator
    del lattice
