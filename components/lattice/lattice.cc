#include "lattice.h"

#include <iostream>
#include <stdexcept>

Lattice::Lattice(std::string const& name) :
    name(name), reference_particle_allocated(false), elements(),
            element_adaptor_map_sptr(new Element_adaptor_map)
{
}

Lattice::Lattice(std::string const& name,
        Element_adaptor_map_sptr const& element_adaptor_map_sptr) :
    name(name), reference_particle_allocated(false), elements(),
            element_adaptor_map_sptr(element_adaptor_map_sptr)
{
}

std::string const&
Lattice::get_name() const
{
    return name;
}

Element_adaptor_map &
Lattice::get_element_adaptor_map()
{
    return *element_adaptor_map_sptr;
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
    Lattice_element_sptr element_sptr(new Lattice_element(element));
    if (element_adaptor_map_sptr->has_adaptor(element.get_type())) {
        element_adaptor_map_sptr->get_adaptor(element.get_type())->set_default_attributes(
                *element_sptr);
    }
    elements.push_back(element_sptr);
}

Lattice_elements &
Lattice::get_elements()
{
    return elements;
}

double
Lattice::get_length() const
{
    double length = 0.0;
    for (Lattice_elements::const_iterator it = elements.begin(); it
            != elements.end(); ++it) {
        length += (*it)->get_length();
    }
    return length;
}

double
Lattice::get_total_angle() const
{
    double angle = 0.0;
    for (Lattice_elements::const_iterator it = elements.begin(); it
            != elements.end(); ++it) {
        angle += (*it)->get_bend_angle();
    }
    return angle;
}

void
Lattice::print() const
{
    std::cout << name << ":\n";
    for (Lattice_elements::const_iterator it = elements.begin(); it
            != elements.end(); ++it) {
        (*it)->print();
    }
}
Lattice::~Lattice()
{
    if (reference_particle_allocated) {
        delete reference_particle_ptr;
    }

}
