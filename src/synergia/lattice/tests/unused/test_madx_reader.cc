#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/lattice/madx_reader.h"
#include "synergia/foundation/physical_constants.h"

const double tolerance = 1.0e-12;

std::string
get_fodo()
{
    std::string fodo("lq=1.0;\n");
    fodo += "ld=2.0;\n";
    fodo += "endmark: marker;\n";
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
get_fodo_with_from()
{
    std::string fodo("lq=1.0;\n");
    fodo += "ld=2.0;\n";
    fodo += "endmark: marker;\n";
    fodo += "f: quadrupole, l=lq;\n";
    fodo += "d: quadrupole, l=lq;\n";
    fodo += "fodo: sequence, l=2*lq+2*ld;\n";
    fodo += "m: marker, at=1.0;\n";
    fodo += "f1: f, at=1;\n";
    fodo += "d1: d, at=3, from=m;\n";
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
    fodo_bodo += "endmark: marker;\n";
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

BOOST_AUTO_TEST_CASE(get_sequence_names)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo_bodo());
    std::vector < std::string > sequence_names(madx_reader.get_sequence_names());
    BOOST_CHECK(sequence_names.size() == 2);
    int fodo_count = 0;
    int bodo_count = 0;
    for (int i = 0; i < 2; ++i) {
        if (sequence_names.at(i) == "fodo") {
            ++fodo_count;
        }
        if (sequence_names.at(i) == "bodo") {
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

BOOST_AUTO_TEST_CASE(get_double_variable)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo());
    const double tolerance = 1.0e-12;
    BOOST_CHECK_CLOSE(madx_reader.get_double_variable("ld"), 2.0, tolerance);
}

BOOST_AUTO_TEST_CASE(get_string_variable)
{
    MadX_reader madx_reader;
    madx_reader.parse("foo='bar';\n" + get_fodo());
    BOOST_CHECK_EQUAL(madx_reader.get_string_variable("foo"), "bar");
}

BOOST_AUTO_TEST_CASE(get_lattice_sptr)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo());
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo"));
//    lattice_sptr->print();
}

BOOST_AUTO_TEST_CASE(get_lattice_sptr2)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo_with_from());
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo"));
//    lattice_sptr->print();

    std::list<Lattice_element_sptr>::const_iterator it = 
        lattice_sptr->get_elements().begin();

    std::advance(it, 3);
    Lattice_element element(**(it));

    const double tolerance = 1.0e-10;
    BOOST_CHECK_CLOSE(element.get_double_attribute("l"), 2.0, tolerance);
}

BOOST_AUTO_TEST_CASE(get_lattice_sptr3)
{
    MadX_reader madx_reader;
    madx_reader.parse_file("lattices/fodo2work.madx");
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo"));
    //lattice_sptr->print();
}

BOOST_AUTO_TEST_CASE(copy_lattice_sptr)
{
    MadX_reader madx_reader;
    madx_reader.parse_file("lattices/fodo2work.madx");
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo"));

    Lattice_sptr copy_lattice_sptr(new Lattice(*lattice_sptr));

    //lattice_sptr->print();
}

BOOST_AUTO_TEST_CASE(get_lattice_sptr_embedded_sequence)
{
    MadX_reader madx_reader;
    madx_reader.parse_file("lattices/cf_oboobo.seq");
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("model"));
    //lattice_sptr->print();

    const double tolerance = 1.0e-10;
    BOOST_CHECK_CLOSE(lattice_sptr->get_length(), 128.0, tolerance);
    BOOST_CHECK_CLOSE(lattice_sptr->get_total_angle(), 6.283185307179588, tolerance);

    Lattice_elements const & elements = lattice_sptr->get_elements();

    int idx = 0;
    for(Lattice_elements::const_iterator it = elements.begin();
            it != elements.end(); ++it, ++idx)
    {
        if (idx == 0)
        {
            BOOST_CHECK( (*it)->get_name().compare("bf") == 0 );
            BOOST_CHECK( (*it)->get_type().compare("sbend") == 0 );
            BOOST_CHECK_CLOSE( (*it)->get_length(), 3, tolerance );
        }
        else if (idx == 1)
        {
            BOOST_CHECK( (*it)->get_name().compare("auto_drift_fodo_0") == 0 );
            BOOST_CHECK( (*it)->get_type().compare("drift") == 0 );
            BOOST_CHECK_CLOSE( (*it)->get_length(), 5, tolerance );
        }
        else if (idx == 2)
        {
            BOOST_CHECK( (*it)->get_name().compare("bd") == 0 );
            BOOST_CHECK( (*it)->get_type().compare("sbend") == 0 );
            BOOST_CHECK_CLOSE( (*it)->get_length(), 3, tolerance );
        }
        else if (idx == 20)
        {
            BOOST_CHECK( (*it)->get_name().compare("bf") == 0 );
            BOOST_CHECK( (*it)->get_type().compare("sbend") == 0 );
            BOOST_CHECK_CLOSE( (*it)->get_length(), 3, tolerance );
        }
        else if (idx == 31)
        {
            BOOST_CHECK( (*it)->get_name().compare("auto_drift_fodo_1") == 0 );
            BOOST_CHECK( (*it)->get_type().compare("drift") == 0 );
            BOOST_CHECK_CLOSE( (*it)->get_length(), 5, tolerance );
        }
    }

}

BOOST_AUTO_TEST_CASE(get_lattice_sptr_no_reference_particle)
{
    MadX_reader madx_reader;
    madx_reader.parse(get_fodo());
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo"));
    BOOST_CHECK(!lattice_sptr->has_reference_particle());
}

BOOST_AUTO_TEST_CASE(get_lattice_sptr_with_reference_particle)
{
    std::string beam_fodo("beam, particle=proton, gamma=1.42;");
    beam_fodo += get_fodo();
    MadX_reader madx_reader;
    madx_reader.parse(beam_fodo);
    Lattice_sptr lattice_sptr(madx_reader.get_lattice_sptr("fodo"));
    BOOST_CHECK(lattice_sptr->has_reference_particle());
    const double tolerance = 1.0e-10;
    BOOST_CHECK_CLOSE(
            lattice_sptr->get_reference_particle().get_mass(),
            pconstants::mp, tolerance);
    BOOST_CHECK_CLOSE(lattice_sptr->get_reference_particle().get_gamma(), 1.42,
            tolerance);
}

BOOST_AUTO_TEST_CASE(get_types)
{
    std::string str(
            "element: quadrupole, x=3.14, name='foo', knl={1.1,2.2,3.3};");
    str += "seq:sequence, l=1.0;\n";
    str += "e1: element, at=0.5;\n";
    str += "endsequence;\n";
    MadX_reader madx_reader;
    madx_reader.parse(str);
//    madx_reader.get_lattice_sptr("seq")->print();

// get the second element
    Lattice_element element(
            **(++madx_reader.get_lattice_sptr("seq")->get_elements().begin()));
    const double tolerance = 1.0e-12;
    BOOST_CHECK_CLOSE(element.get_double_attribute("x"), 3.14, tolerance);
    BOOST_CHECK_EQUAL(element.get_string_attribute("name"), "foo");
    BOOST_CHECK_EQUAL(element.get_vector_attribute("knl").size(), 3);
    BOOST_CHECK_CLOSE(element.get_vector_attribute("knl").at(0), 1.1,
            tolerance);
    BOOST_CHECK_CLOSE(element.get_vector_attribute("knl").at(1), 2.2,
            tolerance);
    BOOST_CHECK_CLOSE(element.get_vector_attribute("knl").at(2), 3.3,
            tolerance);
}

BOOST_AUTO_TEST_CASE(line_length_with_endmark)
{
    std::string str(
            "element: quadrupole, x=3.14, name='foo', knl={1.1,2.2,3.3};");
    str += "endmark: marker;\n";
    str += "seq:sequence, l=1.0;\n";
    str += "e1: element, at=0.5;\n";
    str += "endmark, at=1.0;\n";
    str += "endsequence;\n";
    MadX_reader madx_reader;
    madx_reader.parse(str);
//    madx_reader.get_lattice_sptr("seq")->print();

    const double tolerance = 1.0e-12;
    BOOST_CHECK_CLOSE(madx_reader.get_lattice_sptr("seq")->get_length(), 1.0,
            tolerance);
}

BOOST_AUTO_TEST_CASE(line_length_without_endmark)
{
    std::string str(
            "element: quadrupole, x=3.14, name='foo', knl={1.1,2.2,3.3};");
    str += "seq:sequence, l=1.0;\n";
    str += "e1: element, at=0.5;\n";
    str += "endsequence;\n";
    MadX_reader madx_reader;
    madx_reader.parse(str);
//    madx_reader.get_lattice_sptr("seq")->print();

    const double tolerance = 1.0e-12;
    BOOST_CHECK_CLOSE(madx_reader.get_lattice_sptr("seq")->get_length(), 1.0,
            tolerance);
}

