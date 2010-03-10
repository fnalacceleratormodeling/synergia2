#include <iostream>

#include "operator.h"

Independent_operator::Independent_operator(std::string const& name) :
    Operator(name), slices()
{

}

void
Independent_operator::append_slice(Lattice_element_slice_sptr slice)
{
    slices.push_back(slice);
}

Lattice_element_slices const&
Independent_operator::get_slices() const
{
    return slices;
}

void
Independent_operator::print() const
{
    std::cout << "Independent_operator: " << name << std::endl;
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        (*it)->print();
    }
}

Independent_operator::~Independent_operator()
{

}
