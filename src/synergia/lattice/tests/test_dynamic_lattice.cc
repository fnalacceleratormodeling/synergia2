#include "synergia/utils/catch.hpp"

#include "synergia/lattice/dynamic_lattice.h"
#include "synergia/lattice/madx_reader.h"


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

TEST_CASE("dynamic lattice")
{
    std::string str = R"(

        start = 0.35;

        x = 1.0;

        a: quadrupole, l=0.0, k1=x+1.0;
        b: quadrupole, l=0.2;
        c: quadrupole, l=0.4;
        d: quadrupole, l=0.8;

        seq: sequence, l=3.0;
        a, at=0;
        b, at=0;
        c, at=start+0.1;
        d, at=0.6, from=c;
        endsequence;
    )";
 
    MadX_reader reader;
    reader.parse(str);

    Logger screen(0, LoggerV::DEBUG);

    auto lattice = reader.get_dynamic_lattice("seq");

    CHECK(lattice.get_elements().front().get_double_attribute("k1")
            == Approx(2.0).margin(1e-12));

    lattice.get_lattice_tree().set_variable("x", 3.0);

    CHECK(lattice.get_elements().front().get_double_attribute("k1")
            == Approx(4.0).margin(1e-12));
}
