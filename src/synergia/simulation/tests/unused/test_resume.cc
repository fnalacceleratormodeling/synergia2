#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/resume.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "synergia/bunch/diagnostics_basic.h"
#include "synergia/bunch/diagnostics_full2.h"
#include "synergia/bunch/diagnostics_particles.h"
#include "synergia/bunch/diagnostics_bulk_track.h"
#include "propagator_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
#include "synergia/simulation/lattice_elements_actions.h"

BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(construct, Propagator_fixture)
{
}


// start out just like test_propagator to create a resume context
BOOST_FIXTURE_TEST_CASE(resume, Propagator_fixture)
{
    Bunch_sptr bunch_sptr(new Bunch(l.b.bunch));
    populate_6d(distribution, *bunch_sptr, means, covariances);
    Bunch_simulator bunch_simulator(bunch_sptr);

    Diagnostics_sptr first_step_full2_diag_sptr(
            new Diagnostics_full2("first_full2_per_step.h5"));
    Diagnostics_sptr second_step_full2_diag_sptr(
            new Diagnostics_full2("second_full2_per_step.h5"));

    bunch_simulator.add_per_step(first_step_full2_diag_sptr);
    bunch_simulator.add_per_step(second_step_full2_diag_sptr);

    Diagnostics_sptr first_turn_particles_diag_sptr(
            new Diagnostics_particles("first_particles_per_turn.h5"));
    Diagnostics_sptr second_turn_particles_diag_sptr(
            new Diagnostics_particles("second_particles_per_turn.h5"));
    Diagnostics_sptr turn_tracks_sptr(
                                 new Diagnostics_bulk_track("turn_tracks.h5", 4));

    bunch_simulator.add_per_turn(first_turn_particles_diag_sptr);
    bunch_simulator.add_per_turn(second_turn_particles_diag_sptr);
    bunch_simulator.add_per_turn(turn_tracks_sptr);

    int num_turns = 8;
    int max_turns = 4;
    propagator.propagate(bunch_simulator, num_turns, max_turns);

    // Now read in resume context
    Resume resume;
    resume.propagate(false, 1, false, 1, false, 1);
}


