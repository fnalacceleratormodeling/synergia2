#include "lattice_element_slice.h"
#include "synergia/utils/floating_point.h"

#include <stdexcept>
#include <iostream>
#include <sstream>

const double split_element_tolerance = 1.0e-9;

Lattice_element_slice::Lattice_element_slice(
        Lattice_element const & element) 
    : element(&element)
    , whole(true)
    , left_edge(true)
    , right_edge(true)
    , left(0.0)
    , right(element.get_length())
    , ref_ct(0.0)
{
}

Lattice_element_slice::Lattice_element_slice(
        Lattice_element const & element,
        double l, double r) 
    : element(&element)
    , whole(true)
    , left_edge(true)
    , right_edge(true)
    , left(l)
    , right(r)
    , ref_ct(0.0)
{
    if (left < 0.0)
    {
        throw std::range_error(
                "Lattice_element_slice: left must be >= 0.0");
    }

    if (left < split_element_tolerance) 
    {
        left_edge = true;
        left = 0.0;
    } 
    else 
    {
        left_edge = false;
    }

    double element_length = element.get_length();

    if (right > (element_length + split_element_tolerance)) 
    {
        throw std::range_error(
                "Lattice_element_slice: right must be no greater than "
                "the length of the element");
    }

    if (floating_point_equal(right, element_length, split_element_tolerance)) 
    {
        right_edge = true;
        right = element_length;
    } 
    else 
    {
        right_edge = false;
    }

    whole = (left_edge && right_edge);
}

#if 0
std::string
Lattice_element_slice::as_string() const
{
    std::stringstream sstream;
    if (whole) {
        sstream << "[begin,end] ";
    } else if(left_edge) {
        sstream << "[begin," << right << "] ";
    } else if(right_edge) {
        sstream << "[" << left << ",end] ";
    } else {
        sstream << "[" << left << "," << right << "] ";
    }

    sstream << element_sptr->as_string();
    return sstream.str();
}

void
Lattice_element_slice::print() const
{
    std::cout << as_string() << std::endl;
}
#endif


