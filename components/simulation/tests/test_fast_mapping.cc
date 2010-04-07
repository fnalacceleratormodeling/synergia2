#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/fast_mapping.h"
#include "components/lattice/chef_utils.h"
#include "components/lattice/chef_lattice.h"
#include "components/foundation/reference_particle.h"
#include "components/foundation/physical_constants.h"
#include "components/bunch/populate.h"
#include "utils/multi_array_typedefs.h"
#include "utils/multi_array_print.h"
#include <beamline/beamline_elements.h>
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const int order = 2;
const double coeff = 3.1415;
const double tolerance = 1.0e-12;
const double drift_length = 1.2;
const double quad_length = 0.3;
const double quad_strength = 0.7;
const double mass = constants::mp;
const double real_num = 1.0e11;
const int total_num = 20;
const double total_energy = 125.0;

BOOST_AUTO_TEST_CASE(construct)
{
    Fast_mapping fast_mapping(order);
}

BOOST_AUTO_TEST_CASE(construct2)
{
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);

    // ugh. all this to get a Mapping
    Chef_elements chef_elements;
    ElmPtr drift1(new drift("drift1", drift_length));
    chef_elements.push_back(drift1);
    ElmPtr quad(new quadrupole("quad", quad_length, quad_strength));
    chef_elements.push_back(quad);
    ElmPtr drift2(new drift("drift2", drift_length));
    chef_elements.push_back(drift2);

    JetParticle::createStandardEnvironments(order);
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            reference_particle);
    for (Chef_elements::const_iterator it = chef_elements.begin(); it
            != chef_elements.end(); ++it) {
        (*it)->propagate(jet_particle);
    }

    Fast_mapping fast_mapping(reference_particle, jet_particle.State());
}

BOOST_AUTO_TEST_CASE(add_term)
{
    Fast_mapping fast_mapping(order);

    Fast_mapping_term fast_mapping_term(order);
    int indices[] = { 0, 2, 4 };
    fast_mapping_term.index(0) = indices[0];
    fast_mapping_term.index(1) = indices[1];
    fast_mapping_term.index(2) = indices[2];

    fast_mapping.add_term(0, fast_mapping_term);
}

BOOST_AUTO_TEST_CASE(apply)
{
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);

    // ugh. all this to get a Mapping
    Chef_elements chef_elements;
    ElmPtr drift1(new drift("drift1", drift_length));
    chef_elements.push_back(drift1);
//    ElmPtr quad(new quadrupole("quad", quad_length, quad_strength));
//    chef_elements.push_back(quad);
//    ElmPtr drift2(new drift("drift2", drift_length));
//    chef_elements.push_back(drift2);

    JetParticle::createStandardEnvironments(order);
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            reference_particle);
    for (Chef_elements::const_iterator it = chef_elements.begin(); it
            != chef_elements.end(); ++it) {
        (*it)->propagate(jet_particle);
    }

    Fast_mapping fast_mapping(reference_particle, jet_particle.State());

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

//    multi_array_print(bunch.get_local_particles(), "particles before");
    fast_mapping.apply(bunch);
//    multi_array_print(bunch.get_local_particles(), "particles after");
}
// test_note: We need to check that apply actual produces the correct results.
//            As of this writing, it almost certainly doesn't

BOOST_AUTO_TEST_CASE(write_read_file)
{
    Four_momentum four_momentum(mass, total_energy);
    Reference_particle reference_particle(four_momentum);

    // ugh. all this to get a Mapping
    Chef_elements chef_elements;
    ElmPtr drift1(new drift("drift1", drift_length));
    chef_elements.push_back(drift1);
    ElmPtr quad(new quadrupole("quad", quad_length, quad_strength));
    chef_elements.push_back(quad);
    ElmPtr drift2(new drift("drift2", drift_length));
    chef_elements.push_back(drift2);

    JetParticle::createStandardEnvironments(order);
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            reference_particle);
    for (Chef_elements::const_iterator it = chef_elements.begin(); it
            != chef_elements.end(); ++it) {
        (*it)->propagate(jet_particle);
    }

    Fast_mapping fast_mapping(reference_particle, jet_particle.State());

    fast_mapping.write_to_file("test_fast_mapping.dat");
    Fast_mapping fast_mapping2("test_fast_mapping.dat");
    fast_mapping2.write_to_file("test_fast_mapping2.dat");
}
// test_note: We still need to verify that the Fast_mapping read
//            is the same as the Fast_mapping written.
//            For now, we can only diff the files test_fast_mapping.dat and
//            test_fast_mapping2.dat
