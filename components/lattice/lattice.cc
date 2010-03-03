#include "lattice.h"

#include <iostream>
#include <stdexcept>

Lattice::Lattice(std::string const& name) :
    name(name), reference_particle_allocated(false),
            elements()
{
}

std::string const&
Lattice::get_name() const
{
    return name;
}

void
Lattice::set_reference_particle(Reference_particle const& reference_particle)
{
    reference_particle_ptr = new Reference_particle(reference_particle);
    reference_particle_allocated = true;
}

bool
Lattice::has_reference_particle() const
{
    return reference_particle_allocated;
}

Reference_particle const&
Lattice::get_reference_particle() const
{
    if (!reference_particle_allocated) {
        throw std::runtime_error("Lattice: no reference_particle available");
    }
    return *reference_particle_ptr;
}

void
Lattice::append(Lattice_element const& element)
{
    elements.push_back(element);
}

Lattice_element_list &
Lattice::get_elements()
{
    return elements;
}

double
Lattice::get_length() const
{
    double length = 0.0;
    for (Lattice_element_list::const_iterator it = elements.begin(); it
            != elements.end(); ++it) {
        length += it->get_length();
    }
    return length;
}

double
Lattice::get_total_angle() const
{
    double angle = 0.0;
    for (Lattice_element_list::const_iterator it = elements.begin(); it
            != elements.end(); ++it) {
        angle += it->get_bend_angle();
    }
    return angle;
}

Lattice::~Lattice()
{
    if (reference_particle_allocated) {
        delete reference_particle_ptr;
    }

}
