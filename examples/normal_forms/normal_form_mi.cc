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
#include "synergia/foundation/physical_constants.h"
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

    grid_shape[0] = 16;
    grid_shape[1] = 16;
    grid_shape[2] = 32;
    const int part_per_cell = 1;
    const int num_macro_particles = grid_shape[0] * grid_shape[1]
            * grid_shape[2] * part_per_cell;
    const int seed = 4;
    const double num_real_particles = 1.0e13;
    const int num_steps = 8;
    const int num_turns = 10;

    const int map_order = 1;

    const double trans_emit = 1e-6;

    Lattice_sptr lattice_sptr(new Lattice());
    try {
        xml_load(*lattice_sptr, "mi20-egs-fixed-rf.xml");
    }
    catch (std::runtime_error) {
        std::cerr << "normal_form_mi: failed to find mi20-egs-fixed-rf.xml\n";
    }
#if 0
    lattice_sptr->print();
#endif

    std::cout << "Beam (reference particle) parameters:" << std::endl;
    Reference_particle refpart(lattice_sptr->get_reference_particle());
    std::cout << "    total energy: " << refpart.get_total_energy() << std::endl;
    std::cout << "    momentum:     " << refpart.get_momentum() << std::endl;
    std::cout << "    gamma:        " << refpart.get_gamma() << std::endl;
    std::cout << "    beta:         " << refpart.get_beta() << std::endl;

    std::cout << "Lattice parameters:" << std::endl;
    std::cout << "    length: " << lattice_sptr->get_length() << std::endl;
    std::cout << "    angle: " << lattice_sptr->get_total_angle() << std::endl;
    double freq = 588.0 * refpart.get_beta() * pconstants::c/lattice_sptr->get_length();
    std::cout << "    nominal frequency for harmonic number 588: " << freq <<
      std::endl;

    double total_energy = refpart.get_total_energy();


    Space_charge_3d_open_hockney_sptr space_charge_sptr(
            new Space_charge_3d_open_hockney(Commxx(), grid_shape));

    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    std::cout << "linear normal form checks " <<
      (lattice_simulator.check_linear_normal_form() ? "OK" : "NOT OK") <<
      std::endl;

    boost::shared_ptr<Independent_stepper>  stepper_sptr(new Independent_stepper(lattice_simulator, 1));

    //    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
    //      lattice_simulator, space_charge_sptr, num_steps));
    Propagator propagator(stepper_sptr);

    Commxx comm(MPI_COMM_WORLD);
    Bunch_sptr bunch_sptr(new Bunch(lattice_sptr->get_reference_particle(),
            num_macro_particles, num_real_particles, comm));
    Random_distribution distribution(seed, comm);

    // get actions
    std::vector<double> actions(lattice_simulator.get_stationary_actions(1.44e-3, 3.45e-3, 1.00) );

    std::cout << "Stationary actions: " << actions[0] << ", " << actions[1] << ", " << actions[2] << std::endl;

    if ((actions[0]<0.0) || (actions[1]<0.0) || (actions[2]<0.0)) {
      throw std::runtime_error("unable to satisfy requested moments");
    }

#if 1
    populate_6d_stationary_torus(distribution, *bunch_sptr, actions, lattice_simulator);
#else
    populate_6d_stationary_gaussian(distribution, *bunch_sptr, actions, lattice_simulator);
    #endif

    Diagnostics_sptr per_step_diagnostics(new Diagnostics_basic(
								  "nf_per_step.h5"));

    Diagnostics_sptr per_turn_diagnostics(new Diagnostics_full2(
								  "mi_per_turn.h5"));

    Diagnostics_sptr diag_particles(new Diagnostics_particles("mi_particles_initial.h5", 0, 32768));

    Bunch_simulator bunch_simulator(bunch_sptr);

    bunch_simulator.get_diagnostics_actions().add_per_step(per_step_diagnostics);
    bunch_simulator.get_diagnostics_actions().add_per_turn(per_turn_diagnostics);
    bunch_simulator.get_diagnostics_actions().add_per_turn(diag_particles);

#if 0
    diag_particles_i.set_bunch_sptr(bunch_sptr);
    diag_particles_i.update_and_write();
#endif

    double t0 = MPI_Wtime();
#if 0
    propagator.propagate(*bunch_sptr, num_turns, per_step_diagnostics,
            per_turn_diagnostics, true);
#endif
    propagator.propagate(bunch_simulator, num_turns);

    double t1 = MPI_Wtime();
    if (comm.get_rank() == 0) {
      std::cout << "propagate time = " << (t1-t0) << std::endl;
    }

#if 0
    Diagnostics_particles diag_particles_f(bunch_sptr, "mi_particles_final.h5", 0, 32768);
    diag_particles_f.set_bunch_sptr(bunch_sptr);
    diag_particles_f.update_and_write();
#endif

}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
