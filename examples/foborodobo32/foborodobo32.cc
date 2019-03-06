#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/foundation/reference_particle.h"
#include "synergia/foundation/four_momentum.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/foundation/math_constants.h"
#include "synergia/utils/serialization.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/core_diagnostics.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"
#include "synergia/utils/multi_array_print.h"
#include "synergia/simulation/propagate_actions.h"

#include "ramp_actions.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    Logger logger(0, "foborodobo32.log");

    std::vector<int > grid_shape(3);
    grid_shape[0] = 32;
    grid_shape[1] = 32;
    grid_shape[2] = 256;
    const int part_per_cell = 10;
    const int num_macro_particles = grid_shape[0] * grid_shape[1]
            * grid_shape[2] * part_per_cell;
    const double num_real_particles = 1e13;
    const int num_steps = 8;
    const int num_turns = 4;
    const int map_order = 2;
    const double momentum = 2.0;

    const double stdx = 1.0e-3;
    const double stdy = 1.0e-3;
    //const double stdz = .055707;
    const double stdz = 1.0;

    //const int macro_particles = 10000;
    const int macro_particles = 400;
    const double real_particles = 5.0e10;
    const long int seed = 12345791;

    const int turns = 20;
    const int maxturns=10;
    const int steps = 128;

    Lattice_sptr lattice_sptr(MadX_reader().get_lattice_sptr("model", "foborodobo32.madx"));
    logger << "Read lattice " << lattice_sptr->get_name() << ", length: " << lattice_sptr->get_length() << std::endl;

    Four_momentum four_momentum(pconstants::mp);
    four_momentum.set_momentum(momentum);
    Reference_particle refpart(1, four_momentum);
    double energy = refpart.get_total_energy();
    double beta = refpart.get_beta();
    double gamma = refpart.get_gamma();

    logger << "energy: " << energy << std::endl;
    logger << "momentum: " << momentum << std::endl;
    logger << "gamma: " << gamma << std::endl;
    logger << "beta: " << beta << std::endl;

    lattice_sptr->set_reference_particle(refpart);


    // The RF cavity should be set with a harmonic number = 32

    Lattice_simulator lattice_simulator(lattice_sptr, 1);

    lattice_simulator.register_closed_orbit(false);

    // register closed orbit is supposed to set the cavity frequency
    for (Lattice_elements::const_iterator elemiter = lattice_sptr->get_elements().begin();
         elemiter != lattice_sptr->get_elements().end(); ++elemiter) {
        if ((*elemiter)->get_type() == "rfcavity") {
            logger << "RFCavity: " << (*elemiter)->as_string() << std::endl;
        }
    }

    MArray2d map = lattice_simulator.get_linear_one_turn_map(false);

    double ax,bx,mux;
    double ay,by,muy;
    double as,bs,mus;
    map_to_twiss(map[boost::indices[range(0,2)][range(0,2)]], ax, bx, mux);
    map_to_twiss(map[boost::indices[range(2,4)][range(2,4)]], ay, by, muy);
    map_to_twiss(map[boost::indices[range(4,6)][range(4,6)]], as, bs, mus);

    logger << "lattice parameters ax: " << ax << ", bx: " << bx << ", nu_x: " << mux/(2.0*mconstants::pi) << std::endl;
    logger << "lattice parameters ay: " << ay << ", by: " << by << ", nu_y: " << muy/(2.0*mconstants::pi) << std::endl;
    logger << "lattice parameters as: " << as << ", bs: " << bs << ", nu_s: " << mus/(2.0*mconstants::pi) << std::endl;

    double alpha_c = lattice_simulator.get_momentum_compaction();
    double slip_factor = alpha_c - 1.0/(gamma*gamma);

    logger << "compaction factor: " << alpha_c << std::endl;
    logger << "slip_factor: " << slip_factor << std::endl;

    MArray1d input_means(boost::extents[6]);
    for (int i=0; i<6; ++i) {
        input_means[i] = 0.0;
    }
    MArray2d correlation_matrix(get_correlation_matrix(map, stdx, stdy, stdz, beta));

    Core_diagnostics::print_bunch_parameters(correlation_matrix, beta);

    Commxx_sptr commxx(new Commxx);
    Bunch_sptr bunch_sptr(new Bunch(refpart, macro_particles, real_particles, Commxx_sptr(new Commxx)));

    Random_distribution dist(seed, *commxx);
    populate_6d(dist, *bunch_sptr, input_means, correlation_matrix);

    multi_array_print(bunch_sptr->get_local_particles(), "particles");

    MArray1d bunch_means(Core_diagnostics::calculate_mean(*bunch_sptr));
    MArray1d bunch_stds(Core_diagnostics::calculate_std(*bunch_sptr, bunch_means));

    logger << "Generated bunch statistics:" << std::endl;

    std::string coord_names[] = {"x", "xp", "y", "yp", "cdt", "dpop"};
    for (int i=0; i<6; ++i) {
        logger << coord_names[i] << "\t" << bunch_means[i] << "\t" << bunch_stds[i] << std::endl;
    }

    Bunch_simulator bunch_simulator(bunch_sptr);

    bunch_simulator.add_per_turn(
                Diagnostics_sptr(new Diagnostics_full2("full2.h5")));
    bunch_simulator.add_per_turn(
                Diagnostics_sptr(new Diagnostics_bulk_track("tracks.h5", 100)));
    bunch_simulator.add_per_turn(
                Diagnostics_sptr(new Diagnostics_particles("particles.h5")));


    Split_operator_stepper_sptr stepper_sptr(
                new Split_operator_stepper(lattice_sptr, map_order,
                                           Collective_operator_sptr(new Dummy_collective_operator("foo")),
                                           steps));

    Propagator propagator(stepper_sptr);
    //propagator.set_checkpoint_with_xml(true);
    Ramp_actions ramp_actions(0);
    propagator.propagate(bunch_simulator, ramp_actions, turns, maxturns, 1);

#if 0
    Commxx_sptr commxx_per_host_sptr(new Commxx(true));
    Space_charge_3d_open_hockney_sptr space_charge_sptr(
            new Space_charge_3d_open_hockney(commxx_per_host_sptr, grid_shape));
    space_charge_sptr->set_charge_density_comm(Space_charge_3d_open_hockney::charge_allreduce);
    Lattice_simulator lattice_simulator(lattice_sptr, map_order);
    Split_operator_stepper_sptr stepper_sptr(
            new Split_operator_stepper(lattice_simulator, space_charge_sptr,
                    num_steps));
    Propagator propagator(stepper_sptr);

    Commxx_sptr comm_sptr(new Commxx);
    Bunch_sptr bunch_sptr(
            new Bunch(lattice_sptr->get_reference_particle(),
                    num_macro_particles, num_real_particles, comm_sptr));
    Random_distribution distribution(seed, *comm_sptr);
    MArray1d means;
    xml_load(means, "cxx_means.xml");
    MArray2d covariances;
    xml_load(covariances, "cxx_covariance_matrix.xml");
    populate_6d(distribution, *bunch_sptr, means, covariances);

    Bunch_simulator bunch_simulator(bunch_sptr);
    bunch_simulator.get_diagnostics_actions().add_per_step(
            Diagnostics_sptr(
                    new Diagnostics_basic("cxx_example_per_step.h5")));
    bunch_simulator.get_diagnostics_actions().add_per_turn(
            Diagnostics_sptr(
                    new Diagnostics_full2("cxx_example_per_turn.h5")));
    double t0 = MPI_Wtime();
    const int max_turns = 0;
    const int verbosity = 2;
    Ramp_actions ramp_actions(0);
    propagator.propagate(bunch_simulator, ramp_actions, num_turns, max_turns, verbosity);
    double t1 = MPI_Wtime();
    if (comm_sptr->get_rank() == 0) {
        std::cout << "propagate time = " << (t1 - t0) << std::endl;
    }

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


