#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/lsexpr.h"
#include <iostream>
#include <sstream>

bool
compare_lsexpr_string(Lsexpr const& lsexpr, std::string const& string)
{
    std::stringstream sstream;
    lsexpr.write(sstream);
    bool retval = sstream.str() == string;
    if (!retval) {
        std::cerr << "compare_lsexpr_string():\n\"" <<
                     sstream.str() << "\"\n!=\n\"" <<
                     string << "\"" <<  std::endl;
    }    return retval;
}

bool
compare_lsexpr_lsexpr(Lsexpr const& lsexpr1, Lsexpr const& lsexpr2)
{
    std::stringstream sstream1;
    lsexpr1.write(sstream1);
    std::stringstream sstream2;
    lsexpr2.write(sstream2);
    bool retval = (sstream1.str() == sstream2.str());
    if (!retval) {
        std::cerr << "compare_lsexpr_lsexpr():\n\"" <<
                     sstream1.str() << "\"\n!=\n\"" <<
                     sstream2.str() << "\"" <<  std::endl;
    }    return retval;
}

BOOST_AUTO_TEST_CASE(construct)
{
    Lsexpr lsexpr;
}

BOOST_AUTO_TEST_CASE(write_lsexpr_atom)
{
    Lsexpr lsexpr("atomic_foo");
    BOOST_CHECK(compare_lsexpr_string(lsexpr, "atomic_foo"));
}

BOOST_AUTO_TEST_CASE(write_lsexpr_labeled_atom)
{
    Lsexpr lsexpr("atomic_foo");
    lsexpr.set_label("label");
    BOOST_CHECK(compare_lsexpr_string(lsexpr,
                                      "label: atomic_foo"));
}

BOOST_AUTO_TEST_CASE(write_lsexpr_sequence)
{
    Lsexpr lsexpr;
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.push_back(baz);
    BOOST_CHECK(compare_lsexpr_string(lsexpr,
                                      "{foo, bar, baz}"));
}

BOOST_AUTO_TEST_CASE(write_lsexpr_labeled_sequence)
{
    Lsexpr lsexpr;
    lsexpr.set_label("other_label");
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.push_back(baz);
    BOOST_CHECK(compare_lsexpr_string(lsexpr,
                                      "other_label: \n     {foo, bar, baz}"));
}

BOOST_AUTO_TEST_CASE(write_lsexpr_complex)
{
    Lsexpr lsexpr;
    lsexpr.set_label("complex");
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("bar");
    bar.set_label("second_element");
    lsexpr.push_back(bar);
    Lsexpr weird;
    weird.set_label("weird_elements");
    Lsexpr corge("corge");
    weird.push_back(corge);
    Lsexpr garply("garply");
    weird.push_back(garply);
    lsexpr.push_back(weird);
    BOOST_CHECK(compare_lsexpr_string(lsexpr,
                                      "complex: \n     {foo, second_element: bar, \n     weird_elements: \n         {corge, garply}}"));
}

BOOST_AUTO_TEST_CASE(read_lsexpr_atom)
{
    Lsexpr lsexpr("atomic_foo");

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(read_lsexpr_sequence)
{
    Lsexpr lsexpr;
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.push_back(baz);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(read_lsexpr_labeled_sequence)
{
    Lsexpr lsexpr;
    lsexpr.set_label("other_label");
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.push_back(baz);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(read_lsexpr_complex)
{
    Lsexpr lsexpr;
    lsexpr.set_label("complex");
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("bar");
    bar.set_label("second_element");
    lsexpr.push_back(bar);
    Lsexpr weird;
    weird.set_label("weird_elements");
    Lsexpr corge("corge");
    weird.push_back(corge);
    Lsexpr garply("garply");
    weird.push_back(garply);
    lsexpr.push_back(weird);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(read_quoted_whitespace)
{
    Lsexpr lsexpr;
    Lsexpr foo("foo");
    lsexpr.push_back(foo);
    Lsexpr bar("white space");
    lsexpr.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.push_back(baz);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(int_constructor)
{
    Lsexpr lsexpr(7);
    lsexpr.set_label("int");

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(double_constructor)
{
    Lsexpr lsexpr;
    Lsexpr first(3);
    first.set_label("really int");
    lsexpr.push_back(first);
    Lsexpr second(3.14);
    second.set_label("short");
    lsexpr.push_back(second);
    Lsexpr third(3.141592653689793);
    third.set_label("long");
    lsexpr.push_back(third);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(string_vector_constructor)
{
    std::vector<std::string> string_vector;
    string_vector.push_back("foo bar");
    string_vector.push_back("baz");
    string_vector.push_back("garply");
    Lsexpr lsexpr(string_vector);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(int_vector_constructor)
{
    std::vector<int> int_vector;
    int_vector.push_back(0);
    int_vector.push_back(1);
    int_vector.push_back(1);
    int_vector.push_back(2);
    int_vector.push_back(3);
    int_vector.push_back(5);
    int_vector.push_back(8);
    int_vector.push_back(13);
    Lsexpr lsexpr(int_vector);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(double_vector_constructor)
{
    std::vector<double> double_vector;
    double_vector.push_back(-1.7e32);
    double_vector.push_back(3.141592653589793);
    double_vector.push_back(7.1e129);
    Lsexpr lsexpr(double_vector);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(fake_reference_particle)
{
    Lsexpr refpart;
    refpart.set_label("reference_particle");
    Lsexpr fourmom;
    fourmom.set_label("four_momentum");
    fourmom.push_back(Lsexpr(1.1, "mass"));
    fourmom.push_back(Lsexpr(1.2, "gamma"));
    refpart.push_back(fourmom);
    refpart.push_back(Lsexpr(-1,"charge"));

    std::stringstream sstream;
    refpart.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(refpart, parsed_lsexpr));
}

BOOST_AUTO_TEST_CASE(nested)
{
    Lsexpr lsexpr;
    Lsexpr a;
    Lsexpr a1;
    Lsexpr a11("a11");
    Lsexpr a12("a12");
    a1.push_back(a11);
    a1.push_back(a12);
    Lsexpr a2("a2");
    a.push_back(a1);
    a.push_back(a2);
    lsexpr.push_back(a);
    Lsexpr bar("bar");
    lsexpr.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.push_back(baz);

    std::stringstream sstream;
    lsexpr.write(sstream);
    Lsexpr parsed_lsexpr(sstream);
    BOOST_CHECK(compare_lsexpr_lsexpr(lsexpr, parsed_lsexpr));
}
