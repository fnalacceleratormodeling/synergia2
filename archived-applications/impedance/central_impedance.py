#!/usr/bin/env python
import sys
import synergia
from mpi4py import MPI
import numpy as np
import matplotlib.pyplot as plt
import h5py

def propagate():
    lattice = synergia.lattice.Lattice("foo", synergia.lattice.MadX_adaptor_map())
    dr = synergia.lattice.Lattice_element("drift", "dr")
    dr.set_double_attribute("l", 1.0)
    lattice.append(dr)

    # 400 MeV Booster like
    mp = synergia.foundation.pconstants.mp
    ke = 0.4
    energy = ke+mp
    refpart = synergia.foundation.Reference_particle(1, mp, energy)
    lattice.set_reference_particle(refpart)
    gamma = refpart.get_gamma()
    beta = refpart.get_beta()

    print('energy: ', energy)
    print('gamma: ', gamma)
    print('beta: ', beta)

    orbit_length = 474.202752 # Booster
    harmonic_number = 84
    bucket_length = orbit_length/harmonic_number
    halfbucket = bucket_length/2
    steps = 1
    order = 1
    num_bunches = 1
    verbosity = 20
    turns = 1
    maxturns = 10

    print('bucket length: ', bucket_length)
    npart = 50000
    bunch_charge = 4.5e12/81 # Booster bunch charge

    parent_comm = synergia.utils.Commxx()
    comms = synergia.utils.generate_subcomms(parent_comm, num_bunches)

    bunches = []

    for i in range(num_bunches):
        bunch = synergia.bunch.Bunch(refpart, npart, bunch_charge, comms[i])
        bunch.set_longitudinal_aperture_length(bucket_length)
        bunch.set_bucket_index(0)

        lp = bunch.get_local_particles()
        # Initialize all particles but one at the head of the bunch with a
        # transverse offset

        dx = -0.001
        dy = -0.001
        lp[:, :6] = 0.0
        lp[2:, 0] = dx
        lp[2:, 2] = dy
        # Make triangular distribution between -3/4 and -1/4 of the bucket
        NN = (npart-2)//2
        print('NN: ', NN)
        zmin = -0.5 * halfbucket/beta
        zmax = +0.5 * halfbucket/beta
        zmid = 0.0
        Lhalf = zmid - zmin
        print('zmin: ', zmin)
        print('zmid: ', zmid)
        print('zmax: ', zmax)
        print('Lhalf: ', Lhalf)
        lp[2:NN+2, 4] = zmin + np.sqrt(Lhalf**2 * np.arange(NN)/NN)
        lp[NN+2:2*NN+2, 4] = zmax - np.sqrt(Lhalf**2 * np.arange(NN)/NN)

        # one particle upstream of everything so we don't get messed up
        # by longitudinal binning
        lp[1, 4] = -0.99 * halfbucket/beta # (cT<0 means leading)
  
        # particle at 0 position starts at everything 0

        bunches.append(bunch)

    print('created bunches list, nbunches: ', len(bunches), flush=True)
    bunch_train = synergia.bunch.Bunch_train(bunches, bucket_length)
    print('created bunch_train', flush=True)
    bunch_train_simulator = synergia.simulation.Bunch_train_simulator(bunch_train)
    print('created bunch_train_simulator', flush=True)

    zgrid=1000
    waketype = "XLXTYLYTZpp"
    wave_number = [0, 0, 0]
    full_machine = False
    registered_turns = 1

    operators = []
    imped = synergia.collective.Impedance("Dwake.dat", waketype, zgrid, orbit_length, bucket_length, registered_turns, full_machine, wave_number)
    print('dir(imped): ', dir(imped))
    print('imped name, type: ', imped.get_name(), ',', imped.get_type())
    print(imped)
    operators.append(synergia.simulation.Dummy_collective_operator('foo'))
    operators.append(imped)

    print('created Impedance operator', flush=True)

    print(f'operators ({len(operators)}): ', operators, flush=True)

    stepper = synergia.simulation.Split_operator_stepper(
                                lattice, order, operators, steps)
    print('created stepper', flush=True)

    for i in range(num_bunches):
        bunch_train_simulator.add_per_turn(i, synergia.bunch.Diagnostics_bulk_track(f"central_tracks_{i:02d}.h5",
                                                                           4))
        bunch_train_simulator.add_per_turn(i, synergia.bunch.Diagnostics_full2(f"central_diag_{i:02d}.h5"))

        bunch_train_simulator.add_per_turn(i, synergia.bunch.Diagnostics_particles(f"central_particles_{i:02d}.h5"))

    print('added diagnostics to bunch_train_simulator', flush=True)

    propagator = synergia.simulation.Propagator(stepper)

    propagator.propagate(bunch_train_simulator, turns, maxturns, verbosity)

if __name__ == "__main__":
    propagate()

    h5 = h5py.File('central_tracks_00.h5', 'r')
    tracks = h5.get('track_coords')
    print('particle 0 dpx: ', tracks[1, 0, 1])
    print('particle 0 dpy: ', tracks[1, 0, 3])

    h5.close()

    #try:
    #    main()

    #except Exception as e:
    #    sys.stderr.write(str(e) + '\n')
    #    MPI.COMM_WORLD.Abort(777)
