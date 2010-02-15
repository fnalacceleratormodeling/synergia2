#include "lattice_element.h"

Lattice_element::Lattice_element(std::string const& type,
        std::string const& name) :
    type(type), name(name), double_attributes(), string_attributes()
{
}

std::string const &
Lattice_element::get_type() const
{
    return type;
}

std::string const &
Lattice_element::get_name() const
{
    return name;
}

void
Lattice_element::set_double_attribute(std::string const& name, double value)
{
    double_attributes[name] = value;
}

bool
Lattice_element::has_double_attribute(std::string const& name) const
{
    return (double_attributes.count(name) > 0);
}

double
Lattice_element::get_double_attribute(std::string const& name)
{
    return double_attributes[name];
}
