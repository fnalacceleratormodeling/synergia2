#include "lattice_element_slice.h"
#include "synergia/utils/floating_point.h"

#include <stdexcept>
#include <iostream>

const double split_element_tolerance = 1.0e-9;

Lattice_element_slice::Lattice_element_slice(Lattice_element_sptr lattice_element_sptr) :
    element_sptr(lattice_element_sptr), whole(true), left_edge(true), right_edge(
            true), left(0.0)
{
    right = lattice_element_sptr->get_length();
}

Lattice_element_slice::Lattice_element_slice(Lattice_element_sptr lattice_element_sptr,
        double left, double right) :
    element_sptr(lattice_element_sptr), left(left), right(right)
{
    if (left < 0.0) {
        throw std::range_error("Lattice_element_slice: left must be >= 0.0");
    }
    if (left < split_element_tolerance) {
        left_edge = true;
        left = 0.0;
    } else {
        left_edge = false;
    }

    double element_length = lattice_element_sptr->get_length();
    if (right > (element_length + split_element_tolerance)) {
        throw std::range_error(
                "Lattice_element_slice: right must be no greater than the length of the element");
    }
    if (floating_point_equal(right, element_length, split_element_tolerance)) {
        right_edge = true;
        right = element_length;
    } else {
        right_edge = false;
    }
    if (left_edge && right_edge) {
        whole = true;
    } else {
        whole = false;
    }
}

Lattice_element_slice::Lattice_element_slice()
{
}

bool
Lattice_element_slice::is_whole() const
{
    return whole;
}

bool
Lattice_element_slice::has_left_edge() const
{
    return left_edge;
}

bool
Lattice_element_slice::has_right_edge() const
{
    return right_edge;
}

double
Lattice_element_slice::get_left() const
{
    return left;
}

double
Lattice_element_slice::get_right() const
{
    return right;
}

const Lattice_element &
Lattice_element_slice::get_lattice_element() const
{
    return *element_sptr;
}

Lattice_element &
Lattice_element_slice::get_lattice_element()
{
    return *element_sptr;
}

void
Lattice_element_slice::print() const
{
    if (whole) {
        std::cout << "[begin,end] ";
    } else if(left_edge) {
        std::cout << "[begin," << right << "] ";
    } else if(right_edge) {
        std::cout << "[" << left << ",end] ";
    } else {
        std::cout << "[" << left << "," << right << "] ";
    }

    element_sptr->print();
}

template<class Archive>
    void
    Lattice_element_slice::serialize(Archive & ar, const unsigned int version)
    {
        ar & BOOST_SERIALIZATION_NVP(element_sptr);
        ar & BOOST_SERIALIZATION_NVP(whole);
        ar & BOOST_SERIALIZATION_NVP(left_edge);
        ar & BOOST_SERIALIZATION_NVP(right_edge);
        ar & BOOST_SERIALIZATION_NVP(left);
        ar & BOOST_SERIALIZATION_NVP(right);
    }

template
void
Lattice_element_slice::serialize<boost::archive::binary_oarchive >(
        boost::archive::binary_oarchive & ar, const unsigned int version);

template
void
Lattice_element_slice::serialize<boost::archive::xml_oarchive >(
        boost::archive::xml_oarchive & ar, const unsigned int version);

template
void
Lattice_element_slice::serialize<boost::archive::binary_iarchive >(
        boost::archive::binary_iarchive & ar, const unsigned int version);

template
void
Lattice_element_slice::serialize<boost::archive::xml_iarchive >(
        boost::archive::xml_iarchive & ar, const unsigned int version);
