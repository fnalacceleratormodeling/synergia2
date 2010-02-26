#include "lattice.h"

#include <iostream>
#include <stdexcept>

Lattice::Lattice(std::string const& name) :
    name(name), reference_particle_allocated(false),
            the_elements()
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

double
Lattice::get_length() const
{
    double length = 0.0;
    for (std::list<Lattice_element >::const_iterator it = the_elements.begin(); it
            != the_elements.end(); ++it) {
        length += it->get_length();
    }
    return length;
}

double
Lattice::get_total_angle() const
{
    double angle = 0.0;
    for (std::list<Lattice_element >::const_iterator it = the_elements.begin(); it
            != the_elements.end(); ++it) {
        angle += it->get_bend_angle();
    }
    return angle;
}

void
Lattice::append(Lattice_element & element)
{
    the_elements.push_back(element);
}

std::list<Lattice_element > &
Lattice::elements()
{
    return the_elements;
}

Lattice::~Lattice()
{
    if (!reference_particle_allocated) {
        delete reference_particle_ptr;
    }

}
