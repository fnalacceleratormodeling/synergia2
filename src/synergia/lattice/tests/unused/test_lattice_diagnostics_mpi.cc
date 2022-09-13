#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/lattice_diagnostics.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
const char attribute[] = "test";

struct Fixture
{
    Fixture() :
        lattice_sptr(new Lattice("lattice")), f("quadrupole", "f"),
                o("drift", "o"), d("quadrupole", "d")
    {
        lattice_sptr->append(f);
        lattice_sptr->append(o);
        lattice_sptr->append(d);
        lattice_sptr->append(o);
        BOOST_TEST_MESSAGE("setup fixture");
    }
    ~Fixture()
    {
        BOOST_TEST_MESSAGE("teardown fixture");
    }

    Lattice_sptr lattice_sptr;
    Lattice_element f, o, d;
};

BOOST_FIXTURE_TEST_CASE(construct, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics.h5", attribute);
}

BOOST_FIXTURE_TEST_CASE(update_and_write, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics_mpi.h5", attribute);
    double val = 0.1;
    for(Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it!=lattice_sptr->get_elements().end(); ++it){
        (*it)->set_double_attribute(attribute, val);
        val += 1.0;
    }
    lattice_diagnostics.update_and_write();
}
