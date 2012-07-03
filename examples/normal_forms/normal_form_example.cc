#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/simulation/populate_stationary.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    std::vector<int > grid_shape(3);

    grid_shape[0] = 32;
    grid_shape[1] = 32;
    grid_shape[2] = 128;
    const int part_per_cell = 1;
    const int num_macro_particles = grid_shape[0] * grid_shape[1]
            * grid_shape[2] * part_per_cell;
    const int seed = 4;
    const double num_real_particles = 1.0e13;
    const int num_turns = 4;

    const int map_order = 4;

    Lattice_sptr lattice_sptr(new Lattice());
    try {
        xml_load(*lattice_sptr, "cxx_lattice.xml");
    }
    catch (std::runtime_error) {
        std::cerr << "normal_form_example: failed to find cxx_lattice.xml\n";
        std::cerr << "Run normal_form_example.py to generate cxx_lattice.xml\n";
    }
#if 0
    lattice_sptr->print();
#endif
    Commxx_sptr comm_sptr(new Commxx);
    Space_charge_3d_open_hockney_sptr space_charge_sptr(
            new Space_charge_3d_open_hockney(comm_sptr, grid_shape));

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    std::cout << "linear normal form checks " <<
      (lattice_simulator.check_linear_normal_form() ? "OK" : "NOT OK") <<
      std::endl;

    // Check sample normal form coordinates
    const int nsamp = 2;
    MArray2d samplepoints(boost::extents[nsamp][7]);
    for (int i=0; i<nsamp; ++i) {
      for (int j=0; j<7; ++j) {
	samplepoints[i][j] = 0.0;
      }
    }
    // with this lattice, beta_x=33.379, beta_y = 5.867, beta_z =89.31
    // full bucket length = 5.0
    // get particle with maximum extent in x and dp/p
    samplepoints[0][0] = 1.0e-3;
    samplepoints[0][5] = 2.4/89.31/10.0;
    // particle with maximum extent in y
    samplepoints[1][2] = 1.0e-3/2.38; // sqrt(bx/by)

    lattice_simulator.convert_human_to_normal(samplepoints);
    for (int i=0; i<nsamp; ++i) {
      std::cout << "normal form action, angles for particle: " << i << std::endl;
      for (int j=0; j<3; ++j) {
	std::complex<double>z(samplepoints[i][2*j], samplepoints[i][2*j+1]);
	std::cout << "    |a(" << j << ")|: " << std::abs(z) << " ang(a(" << j << ")): " << std::arg(z) << std::endl;
      }
    }

    boost::shared_ptr<Independent_stepper>  stepper_sptr(new Independent_stepper(lattice_simulator, 1));

    //    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
    //      lattice_simulator, space_charge_sptr, num_steps));
    Propagator propagator(stepper_sptr);

    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(),
            num_macro_particles, num_real_particles, comm_sptr));
    Random_distribution distribution(seed, *comm_sptr);

    // get actions

    std::vector<double> actions(lattice_simulator.get_stationary_actions(3e-3, 1e-3, 5e-2) );

    std::cout << "Stationary actions: " << actions[0] << ", " << actions[1] << ", " << actions[2] << std::endl;

    if ((actions[0]<0.0) || (actions[1]<0.0) || (actions[2]<0.0)) {
      throw std::runtime_error("unable to satisfy requested moments");
    }

    std::cout << "before populate" << std::endl;

#if 0
    populate_6d_stationary_torus(distribution, *bunch_sptr, actions, lattice_simulator);
#else
    populate_6d_stationary_gaussian(distribution, *bunch_sptr, actions, lattice_simulator);
    #endif
    std::cout << "after populate" << std::endl;

    Diagnostics_sptr  per_step_diagnostics(new Diagnostics_basic(
            "nf_per_step.h5"));

    Diagnostics_sptr per_turn_diagnostics(new Diagnostics_full2(
							"nf_per_turn.h5"));

    Diagnostics_particles diag_particles_i("nf_particles_initial.h5", 0, 32768);
    diag_particles_i.set_bunch_sptr(bunch_sptr);
    diag_particles_i.update_and_write();

    Bunch_simulator bunch_simulator(bunch_sptr);
    bunch_simulator.add_per_step(per_step_diagnostics);
    bunch_simulator.add_per_turn(per_turn_diagnostics);
    double t0 = MPI_Wtime();
    propagator.propagate(bunch_simulator, num_turns);
    double t1 = MPI_Wtime();
    if (comm_sptr->get_rank() == 0) {
      std::cout << "propagate time = " << (t1-t0) << std::endl;
    }
    Diagnostics_particles diag_particles_f("nf_particles_final.h5", 0, 32768);
    diag_particles_f.set_bunch_sptr(bunch_sptr);
    diag_particles_f.update_and_write();

}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
