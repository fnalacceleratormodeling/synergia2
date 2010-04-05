#include <iostream>

#include "operator.h"

Operator::Operator(std::string const& name) :
    name(name), slices()
{
}

std::string const&
Operator::get_name() const
{
    return name;
}

Lattice_element_slices &
Operator::get_slices()
{
    return slices;
}

void
Operator::apply(Bunch & bunch, Chef_lattice & chef_lattice)
{

}

void
Operator::print() const
{
    std::cout << "Operator " << name << std::endl;
}

Operator::~Operator()
{
}

Collective_operator::Collective_operator(std::string const& name) :
    Operator(name)
{
}

void
Collective_operator::print() const
{
    std::cout << "Collective_operator: " << name << std::endl;
}

Collective_operator::~Collective_operator()
{
}

void
Independent_operator::update_operations(Chef_lattice & chef_lattice)
{
    operations.clear();

    // Group slices of equal operation_type and pass to operation_extractor
    // to get operations.
    Lattice_element_slices group;
    std::string operation_type(""), last_operation_type("");
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        if ((*it)->get_lattice_element().has_string_attribute("operation_type")) {
            operation_type = (*it)->get_lattice_element().get_string_attribute(
                    "operation_type");
        } else {
            operation_type = "default";
        }
        if ((operation_type != last_operation_type)) {
            if (!group.empty()) {
                Independent_operations group_operations =
                        params_ptr->get_extractor(operation_type)->extract(
                                group, chef_lattice);
                for (Independent_operations::const_iterator group_it =
                        group_operations.begin(); group_it
                        != group_operations.end(); ++group_it) {
                    operations.push_back(*group_it);
                }
            }
        }
    }

    have_operations = true;
}

bool
Independent_operator::need_update()
{
    // jfa: this is a placeholder to be replaced when the update mechanism is in place
    return !have_operations;
}

Independent_operator::Independent_operator(std::string const& name,
        Independent_params * ind_params) :
    Operator(name), have_operations(false), params_ptr(ind_params)
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
Independent_operator::apply(Bunch & bunch, Chef_lattice & chef_lattice)
{
    if (need_update()) {
        update_operations(chef_lattice);
    }
    for (Independent_operations::iterator it = operations.begin(); it
            != operations.end(); ++it) {
        (*it)->apply(bunch);
    }
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
