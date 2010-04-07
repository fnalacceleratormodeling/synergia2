#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/fast_mapping.h"
#include "components/foundation/physical_constants.h"
#include "components/bunch/bunch.h"
#include "utils/boost_test_mpi_fixture.h"
#include "components/bunch/populate.h"
#include "utils/multi_array_typedefs.h"
#include "utils/multi_array_print.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double mass = constants::mp;
const double real_num = 1.0e11;
const int total_num = 20;
const double total_energy = 125.0;
const double tolerance = 1.0e-12;
const double drift_length = 1.2;
const double quad_length = 0.3;
const double quad_strength = 0.7;

BOOST_AUTO_TEST_CASE(construct)
{
    Chef_elements chef_elements;
    ElmPtr drift1(new drift("drift1", drift_length));
    chef_elements.push_back(drift1);
    ElmPtr quad(new quadrupole("quad", quad_length, quad_strength));
    chef_elements.push_back(quad);
    ElmPtr drift2(new drift("drift2", drift_length));
    chef_elements.push_back(drift2);

    Chef_propagator chef_propagator(chef_elements);
}

BOOST_AUTO_TEST_CASE(apply)
{
    Chef_elements chef_elements;
    ElmPtr drift1(new drift("drift1", drift_length));
    chef_elements.push_back(drift1);
    ElmPtr quad(new quadrupole("quad", quad_length, quad_strength));
    chef_elements.push_back(quad);
    ElmPtr drift2(new drift("drift2", drift_length));
    chef_elements.push_back(drift2);

    Chef_propagator chef_propagator(chef_elements);

    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);

    Commxx comm(MPI_COMM_WORLD);
    Bunch bunch(reference_particle, constants::proton_charge, total_num,
            real_num, comm);

    MArray2d covariances(boost::extents[6][6]);
    MArray1d means(boost::extents[6]);
    for (int i = 0; i < 6; ++i) {
        means[i] = 0.0;
        for (int j = i; j < 6; ++j) {
            covariances[i][j] = covariances[j][i] = (i + 1) * (j + 1);
        }
    }
    Random_distribution distribution(0, comm);
    populate_6d(distribution, bunch, means, covariances);

//    multi_array_print(bunch.get_local_particles(),"particles before");
    chef_propagator.apply(bunch);
//    multi_array_print(bunch.get_local_particles(),"particles after");
}

// test_note: the apply test just verifies that the method doesn't crash.
//            At least one quantitative test would be desirable.
