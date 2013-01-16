#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/madx_reader.h"

const double tolerance = 1.0e-12;

std::string
get_fodo()
{
    std::string fodo("lq=1.0;\n");
    fodo += "ld=2.0;\n";
    fodo += "f: quadrupole, l=lq;\n";
    fodo += "d: quadrupole, l=lq;\n";
    fodo += "fodo: sequence, l=2*lq+2*ld;\n";
    fodo += "f1: f, at=0;\n";
    fodo += "d1: d, at=3;\n";
    fodo += "endmark, at=2*lq+2*ld;\n";
    fodo += "endsequence;\n";
    return fodo;
}

std::string
get_fodo_bodo()
{
    std::string fodo_bodo(get_fodo());
    fodo_bodo += "lb:=3.0;\n";
    fodo_bodo += "lq2:=3.0;\n";
    fodo_bodo += "ld2:=2.0;\n";
    fodo_bodo += "b: sbend, l=lb;\n";
    fodo_bodo += "bodo: sequence, l=lb+lq2+2*ld2;\n";
    fodo_bodo += "b1: b, at=0;\n";
    fodo_bodo += "d1: d, at=lb+ld2, l=lq2;\n";
    fodo_bodo += "endmark, at=lb+lq2+2*ld;\n";
    fodo_bodo += "endsequence;\n";
    return fodo_bodo;
}

BOOST_AUTO_TEST_CASE(construct)
{
    MadX_reader madx_reader;
}

BOOST_AUTO_TEST_CASE(parse)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo());
}

BOOST_AUTO_TEST_CASE(parse_file)
{
    MadX_reader madx_reader;
    madx_reader.parse_file("lattices/fodo.madx");
}

BOOST_AUTO_TEST_CASE(get_line_names)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo_bodo());
    std::vector < std::string > line_names(madx_reader.get_line_names());
    BOOST_CHECK(line_names.size() == 2);
    int fodo_count = 0;
    int bodo_count = 0;
    for (int i = 0; i < 2; ++i) {
        if (line_names.at(i) == "fodo") {
            ++fodo_count;
        }
        if (line_names.at(i) == "bodo") {
            ++bodo_count;
        }
    }
    BOOST_CHECK(fodo_count == 1);
    BOOST_CHECK(bodo_count == 1);
}

BOOST_AUTO_TEST_CASE(get_line_names_bad)
{
    MadX_reader madx_reader;
    bool caught(false);
    try {
        std::vector < std::string > line_names(madx_reader.get_line_names());
    }
    catch (std::runtime_error) {
        caught = true;
    }
    BOOST_CHECK(caught);
}

BOOST_AUTO_TEST_CASE(get_lattice)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo());
    Lattice lattice(madx_reader.get_lattice("fodo"));
    lattice.print();
}

BOOST_AUTO_TEST_CASE(get_lattice3)
{
    MadX_reader madx_reader;
    madx_reader.parse_file("lattices/fodo2work.madx");
    Lattice lattice(madx_reader.get_lattice("fodo"));
    lattice.print();
}

