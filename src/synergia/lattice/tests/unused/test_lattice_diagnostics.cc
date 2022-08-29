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

BOOST_FIXTURE_TEST_CASE(set_get_default_value, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics.h5", attribute);
    double value = 7.89;
    lattice_diagnostics.set_default_value(value);
    BOOST_CHECK_CLOSE(lattice_diagnostics.get_default_value(), value,
            tolerance);
}

BOOST_FIXTURE_TEST_CASE(set_get_reduce, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics.h5", attribute);
    lattice_diagnostics.set_reduce(true);
    BOOST_CHECK(lattice_diagnostics.get_reduce());

    lattice_diagnostics.set_reduce(false);
    BOOST_CHECK(!lattice_diagnostics.get_reduce());
}

BOOST_FIXTURE_TEST_CASE(set_get_reduce_op, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics.h5", attribute);
    MPI_Op op(MPI_SUM);
    lattice_diagnostics.set_reduce_op(op);
    BOOST_CHECK_EQUAL(lattice_diagnostics.get_reduce_op(), op);

    op = MPI_MAX;
    lattice_diagnostics.set_reduce_op(op);
    BOOST_CHECK_EQUAL(lattice_diagnostics.get_reduce_op(), op);
}

BOOST_FIXTURE_TEST_CASE(is_serial, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics.h5", attribute);
    BOOST_CHECK(lattice_diagnostics.is_serial());
}

BOOST_FIXTURE_TEST_CASE(update_and_write, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics_real.h5", attribute);
    double val = 0.1;
    for(Lattice_elements::iterator it = lattice_sptr->get_elements().begin();
            it!=lattice_sptr->get_elements().end(); ++it){
        (*it)->set_double_attribute(attribute, val);
        val += 1.0;
    }
    lattice_diagnostics.update_and_write();
}

BOOST_FIXTURE_TEST_CASE(update_and_write_default, Fixture)
{
    Lattice_diagnostics lattice_diagnostics(lattice_sptr,
            "lattice_diagnostics_default.h5", attribute);
    lattice_diagnostics.update_and_write();
}
