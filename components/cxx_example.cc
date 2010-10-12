#include <iostream>

// n.b. the use of relative pathnames here is not for general consumption
#include "lattice/lattice.h"
#include "../utils/xml_serialization.h"
#include "simulation/operator.h"
#include "simulation/lattice_simulator.h"
#include "simulation/stepper.h"
#include "simulation/propagator.h"
#include "bunch/bunch.h"
#include "foundation/distribution.h"
#include "bunch/populate.h"
#include "bunch/diagnostics_writer.h"

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);

    const int num_macro_particles = 32000;
    const int seed = 4;
    //    grid = [16, 16, 16]
    const double num_real_particles = 1e12;
    const int num_steps = 8;
    const int num_turns = 4;
    const int map_order = 2;
    const double emit = 1e-6;
    const double stdz = 0.01;
    const double dpop = 1e-4;

    Lattice_sptr lattice_sptr(new Lattice());
    xml_load(*lattice_sptr, "cxx_lattice.xml");
    Collective_operator_sptr space_charge_sptr(new Collective_operator(
            "space_charge"));
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge_sptr, num_steps));
    Propagator propagator(stepper_sptr);

    Commxx comm(MPI_COMM_WORLD);
    Bunch bunch(lattice_sptr->get_reference_particle(), num_macro_particles,
            num_real_particles, comm);
    Random_distribution distribution(seed, comm);
    MArray1d means;
    xml_load(means, "cxx_means.xml");
    MArray2d covariances;
    xml_load(covariances, "cxx_covariance_matrix.xml");
    populate_6d(distribution, bunch, means, covariances);

    Diagnostics_sptr diagnostics_sptr(new Diagnostics);
    Diagnostics_writer per_step_diagnostics("cxx_example_per_step.h5",
            diagnostics_sptr);
    Diagnostics_full2_sptr diagnostics_full2_sptr(new Diagnostics_full2);
    Diagnostics_writer per_turn_diagnostics("cxx_example_per_turn.h5",
            diagnostics_full2_sptr);
    propagator.propagate(bunch, num_turns, per_step_diagnostics,
            per_turn_diagnostics);

    MPI_Finalize();
    return 0;
}
