#define BOOST_TEST_MAIN
#include <boost/test/unit_test.hpp>
#include "synergia/utils/lsexpr.h"
#include <iostream>

BOOST_AUTO_TEST_CASE(construct)
{
    Lsexpr lsexpr;
}

BOOST_AUTO_TEST_CASE(write_lsexpr_atom)
{
    Lsexpr lsexpr("atomic_foo");
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(write_lsexpr_labeled_atom)
{
    Lsexpr lsexpr("atomic_foo");
    lsexpr.set_label("label");
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(write_lsexpr_sequence)
{
    Lsexpr lsexpr;
    Lsexpr foo("foo");
    lsexpr.sequence.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.sequence.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.sequence.push_back(baz);
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(write_lsexpr_labeled_sequence)
{
    Lsexpr lsexpr;
    lsexpr.set_label("other_label");
    Lsexpr foo("foo");
    lsexpr.sequence.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.sequence.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.sequence.push_back(baz);
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(write_lsexpr_complex)
{
    Lsexpr lsexpr;
    lsexpr.set_label("complex");
    Lsexpr foo("foo");
    lsexpr.sequence.push_back(foo);
    Lsexpr bar("bar");
    bar.set_label("second_element");
    lsexpr.sequence.push_back(bar);
    Lsexpr weird;
    weird.set_label("weird_elements");
    Lsexpr corge("corge");
    weird.sequence.push_back(corge);
    Lsexpr garply("garply");
    weird.sequence.push_back(garply);
    lsexpr.sequence.push_back(weird);
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(read_lsexpr_atom)
{
    Lsexpr lsexpr("atomic_foo");
    std::cout << "original:\n";
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
    std::stringstream sstream;
    write_lsexpr(lsexpr, sstream);
    Lsexpr parsed_lsexpr(read_lsexpr(sstream));
    std::cout << "regurgitated:\n";
    write_lsexpr(parsed_lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(read_lsexpr_sequence)
{
    Lsexpr lsexpr;
    Lsexpr foo("foo");
    lsexpr.sequence.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.sequence.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.sequence.push_back(baz);
    std::cout << "original:\n";
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
    std::stringstream sstream;
    write_lsexpr(lsexpr, sstream);
    Lsexpr parsed_lsexpr(read_lsexpr(sstream));
    std::cout << "regurgitated:\n";
    write_lsexpr(parsed_lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(read_lsexpr_labeled_sequence)
{
    Lsexpr lsexpr;
    lsexpr.set_label("other_label");
    Lsexpr foo("foo");
    lsexpr.sequence.push_back(foo);
    Lsexpr bar("bar");
    lsexpr.sequence.push_back(bar);
    Lsexpr baz("baz");
    lsexpr.sequence.push_back(baz);
    std::cout << "original:\n";
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
    std::stringstream sstream;
    write_lsexpr(lsexpr, sstream);
    Lsexpr parsed_lsexpr(read_lsexpr(sstream));
    std::cout << "regurgitated:\n";
    write_lsexpr(parsed_lsexpr, std::cout);
    std::cout << std::endl;
}

BOOST_AUTO_TEST_CASE(read_lsexpr_complex)
{
    Lsexpr lsexpr;
    lsexpr.set_label("complex");
    Lsexpr foo("foo");
    lsexpr.sequence.push_back(foo);
    Lsexpr bar("bar");
    bar.set_label("second_element");
    lsexpr.sequence.push_back(bar);
    Lsexpr weird;
    weird.set_label("weird_elements");
    Lsexpr corge("corge");
    weird.sequence.push_back(corge);
    Lsexpr garply("garply");
    weird.sequence.push_back(garply);
    lsexpr.sequence.push_back(weird);
    std::cout << "original:\n";
    write_lsexpr(lsexpr, std::cout);
    std::cout << std::endl;
    std::stringstream sstream;
    write_lsexpr(lsexpr, sstream);
    Lsexpr parsed_lsexpr(read_lsexpr(sstream));
    std::cout << "regurgitated:\n";
    write_lsexpr(parsed_lsexpr, std::cout);
    std::cout << std::endl;
}
