#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "components/simulation/fast_mapping.h"
#include "chef_elements_fixture.h"
#include "bunch_fixture.h"
#include "fast_mapping_term_fixture.h"
#include "utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture)
;

const double tolerance = 1.0e-12;

struct Mapping_fixture
{
    Mapping_fixture()
    {
        BOOST_TEST_MESSAGE("setup Mapping fixture");
        JetParticle::createStandardEnvironments(order);
        JetParticle jet_particle = reference_particle_to_chef_jet_particle(
                b.reference_particle);
        for (Chef_elements::const_iterator it = c.chef_elements.begin(); it
                != c.chef_elements.end(); ++it) {
            (*it)->propagate(jet_particle);
        }
        mapping = jet_particle.State();
    }
    ~Mapping_fixture()
    {
        BOOST_TEST_MESSAGE("teardown Mapping fixture");
    }
    Bunch_fixture b;
    Chef_elements_fixture c;
    Mapping mapping;
};

BOOST_AUTO_TEST_CASE(construct)
{
    Fast_mapping fast_mapping(order);
}

BOOST_FIXTURE_TEST_CASE(construct2, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping);
}

BOOST_FIXTURE_TEST_CASE(add_term, Fast_mapping_term_fixture)
{
    Fast_mapping fast_mapping(order);
    fast_mapping.add_term(0, fast_mapping_term);
}

BOOST_FIXTURE_TEST_CASE(apply, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping);
    //    multi_array_print(bunch.get_local_particles(), "particles before");
    fast_mapping.apply(b.bunch);
    //    multi_array_print(bunch.get_local_particles(), "particles after");
}
// test_note: We need to check that apply actual produces the correct results.
//            As of this writing, it almost certainly doesn't

BOOST_FIXTURE_TEST_CASE(write_read_file, Mapping_fixture)
{
    Fast_mapping fast_mapping(b.reference_particle, mapping);

    fast_mapping.write_to_file("test_fast_mapping.dat");
    Fast_mapping fast_mapping2("test_fast_mapping.dat");
    fast_mapping2.write_to_file("test_fast_mapping2.dat");
}
// test_note: We still need to verify that the Fast_mapping read
//            is the same as the Fast_mapping written.
//            For now, we can only diff the files test_fast_mapping.dat and
//            test_fast_mapping2.dat
