#include <iostream>
#include <stdexcept>

#include "synergia/utils/synergia_omp.h"

#include "synergia/foundation/distribution.h"
#include "synergia/utils/serialization.h"

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"

#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/independent_stepper_elements.h"
#include "synergia/simulation/propagator.h"

#include "synergia/bunch/bunch.h"
#include "synergia/bunch/populate.h"

#include "benchmark_options.h"

#include <sched.h>

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
    const int num_steps = 1;
    const int num_turns = 1;
    const int map_order = 2;

    MadX_reader reader;
    Lattice_sptr lattice = reader.get_lattice_sptr("fodo", "fodo.madx");

    //lattice->set_all_string_attribute("extractor_type","chef_propagate");
    lattice->set_all_string_attribute("extractor_type","libff");

    Commxx_sptr comm_sptr(new Commxx);

    Bunch_sptr bunch(
            new Bunch(lattice->get_reference_particle(),
                    num_macro_particles, num_real_particles, comm_sptr));

    std::cout << "bunch local particle = " << bunch->get_local_num() << "\n";

    Random_distribution distribution(seed, *comm_sptr);

    MArray1d means;
    xml_load(means, "cxx_means.xml");

    MArray2d covariances;
    xml_load(covariances, "cxx_covariance_matrix.xml");

    populate_6d(distribution, *bunch, means, covariances);

    Bunch_simulator bunch_simulator(bunch);

    Independent_stepper_elements_sptr stepper(
            new Independent_stepper_elements(lattice, map_order, num_steps));

    Propagator propagator(stepper);
    propagator.set_final_checkpoint(false);

    double t0 = MPI_Wtime();
    const int max_turns = 0;
    propagator.propagate(bunch_simulator, num_turns, max_turns, opts.verbosity);
    double t1 = MPI_Wtime();
    if (comm_sptr->get_rank() == 0) {
        std::cout << "propagate time = " << (t1 - t0) << std::endl;
    }
}

int
main(int argc, char **argv)
{
    int thread_level;
    MPI_Init_thread(&argc, &argv, MPI_THREAD_FUNNELED, &thread_level);

    Benchmark_options opts(argc, argv);
    std::cout << "OpenMP num threads = " << opts.ompthreads << "\n";
    omp_set_num_threads(opts.ompthreads);

    int cpu = sched_getcpu();
    std::cout << "cpuid = " << cpu << "\n";

    run(opts);
    MPI_Finalize();
    return 0;
}
