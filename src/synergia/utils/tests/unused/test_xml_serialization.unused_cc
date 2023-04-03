#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/serialization.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/utils/boost_test_mpi_fixture.h"
BOOST_GLOBAL_FIXTURE(MPI_fixture);

const double tolerance = 1.0e-12;
class Foo
{
private:
    int i;
    double x;

public:
    Foo()
    {
    }
    Foo(int i, double x)
    {
        this->i = i;
        this ->x = x;
    }
    bool
    operator==(Foo const& other)
    {
        return ((this->i == other.i) && (this->x == other.x));
    }
    template<class Archive>
        void
        serialize(Archive & ar, const unsigned int version)
        {
            ar & BOOST_SERIALIZATION_NVP(i) & BOOST_SERIALIZATION_NVP(x);
        }
};

BOOST_AUTO_TEST_CASE(load_save)
{
    Foo a(1, 2.2);
    xml_save<Foo > (a, "a.xml");
    Foo b;
    xml_load<Foo > (b, "a.xml");
    BOOST_CHECK(a==b);
}

