#define BOOST_TEST_MAIN
#include <boost/lexical_cast.hpp>
#include <boost/test/unit_test.hpp>
#include "synergia/simulation/propagator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/split_operator_stepper.h"
#include "synergia/foundation/physical_constants.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/bunch/bunch.h"
#include "lattice_fixture.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
#include "synergia/simulation/diagnostics_actions.h"
#include "synergia/bunch/diagnostics_particles.h"


BOOST_GLOBAL_FIXTURE(MPI_fixture)

const double tolerance = 1.0e-12;

BOOST_FIXTURE_TEST_CASE(propagate, Lattice_fixture)
 {
   Dummy_collective_operator_sptr space_charge(new Dummy_collective_operator(
           "space_charge"));

   Lattice_simulator lattice_simulator(lattice_sptr, 2);

    Split_operator_stepper_sptr stepper_sptr(new Split_operator_stepper(
            lattice_simulator, space_charge, 4));
    Propagator propagator(stepper_sptr);

    Reference_particle reference_particle=b.bunch.get_reference_particle() ;
    /// bunches with different comms has to be generated
    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] =0.;// i * 0.0000072;
        for (int j = i; j < 6; ++j) {
            covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1)*0.00000000001;
        }
         covariances[i][i] *= 10.0; // this makes for a positive-definite matrix
    }

    for (int i = 0; i < 6; ++i) {
         covariances[1][i] *=0.01;
         covariances[i][1] *=0.01;
         covariances[3][i] *=0.01;
         covariances[i][3] *=0.01;
         covariances[5][i] *=0.01;
         covariances[i][5] *=0.01;
    }
 }

