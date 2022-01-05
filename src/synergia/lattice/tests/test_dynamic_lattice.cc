#include "synergia/utils/catch.hpp"

#include "synergia/lattice/dynamic_lattice.h"


TEST_CASE("construct")
{
    REQUIRE_NOTHROW(Dynamic_lattice());
}

TEST_CASE("print")
{
    Dynamic_lattice dl;
    dl.set_variable("x", "1");
    dl.set_variable("a", "x+3");
    //dl.print();
}
