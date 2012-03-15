#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"

#include "benchmark_options.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run(Benchmark_options const& opts)
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = opts.gridx;
    grid_shape[1] = opts.gridy;
    grid_shape[2] = opts.gridz;
    const int part_per_cell = opts.partpercell;
    const int num_macro_particles = grid_shape[0] * grid_shape[1]
            * grid_shape[2] * part_per_cell;
    const int seed = 4;
    const double num_real_particles = 1e13;
    const int num_steps = 8;
    const int num_turns = 4;
    const int map_order = 2;
    const double emit = 1e-6;
    const double stdz = 0.01;
    const double dpop = 1e-4;

    Lattice_sptr lattice_sptr(new Lattice());
    try {
        xml_load(*lattice_sptr, "cxx_lattice.xml");
    }
    catch (std::runtime_error) {
        std::cerr << "cxx_example: failed to find cxx_lattice.xml\n";
        std::cerr << "Run cxx_example.py to generate cxx_lattice.xml\n";
        exit(1);
    }
    Space_charge_3d_open_hockney_sptr space_charge_sptr(
            new Space_charge_3d_open_hockney(Commxx(), grid_shape));
    if (opts.autotune) {
        space_charge_sptr->auto_tune_comm(true);
    } else {
        if (opts.chargecomm > 0) {
            space_charge_sptr->set_charge_density_comm(
                    Space_charge_3d_open_hockney::Charge_density_comm(
                            opts.chargecomm));
        }
        if (opts.efieldcomm > 0) {
            space_charge_sptr->set_e_field_comm(
                    Space_charge_3d_open_hockney::E_field_comm(opts.efieldcomm));
        }
    }
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge_sptr, num_steps));
    Propagator propagator(stepper_sptr);

    Commxx comm(MPI_COMM_WORLD);
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(),
            num_macro_particles, num_real_particles, comm));
    Random_distribution distribution(seed, comm);
    MArray1d means;
    xml_load(means, "cxx_means.xml");
    MArray2d covariances;
    xml_load(covariances, "cxx_covariance_matrix.xml");
    populate_6d(distribution, *bunch_sptr, means, covariances);
    if (opts.sortperiod > 0) {
        bunch_sptr->sort(Bunch::z);
    }
    bunch_sptr->set_sort_period(opts.sortperiod);

    Bunch_simulator bunch_simulator(bunch_sptr);
    bunch_simulator.get_diagnostics_actions().add_per_step(
            Diagnostics_sptr(
                    new Diagnostics_basic("cxx_example_per_step.h5")));
    bunch_simulator.get_diagnostics_actions().add_per_turn(
            Diagnostics_sptr(
                    new Diagnostics_full2("cxx_example_per_turn.h5")));
    double t0 = MPI_Wtime();
    const int max_turns = 0;
    propagator.propagate(bunch_simulator, num_turns, max_turns, opts.verbosity);
    double t1 = MPI_Wtime();
    if (comm.get_rank() == 0) {
      std::cout << "propagate time = " << (t1-t0) << std::endl;
    }
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    Benchmark_options opts(argc, argv);
    run(opts);
    MPI_Finalize();
    return 0;
}
