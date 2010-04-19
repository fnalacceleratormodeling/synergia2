#include <iostream>

#include "operator.h"

Operator::Operator(std::string const& name, std::string const& type) :
    name(name), type(type)
{
}

std::string const&
Operator::get_name() const
{
    return name;
}

std::string const&
Operator::get_type() const
{
    return type;
}

void
Operator::print() const
{
    std::cout << type << " operator: " << name << std::endl;
}

Operator::~Operator()
{
}

Collective_operator::Collective_operator(std::string const& name) :
    Operator(name, "collective")
{
}

void
Collective_operator::apply(Bunch & bunch, Operators & step_operators)
{
    std::cout << "stub: Collective_operator::apply\n";
}

Collective_operator::~Collective_operator()
{
}

void
Independent_operator::update_operations(
        Reference_particle const& reference_particle)
{
    operations.clear();

    // Group slices of equal extractor_type and pass to operation_extractor
    // to get operations.
    Lattice_element_slices group;
    std::string extractor_type(""), last_extractor_type("");
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        if ((*it)->get_lattice_element().has_string_attribute("extractor_type")) {
            extractor_type = (*it)->get_lattice_element().get_string_attribute(
                    "extractor_type");
        } else {
            extractor_type = "default";
        }
        if ((extractor_type != last_extractor_type) && (!group.empty())) {
            Independent_operations
                    group_operations =
                            operation_extractor_map_sptr->get_extractor(
                                    extractor_type)->extract(
                                    reference_particle, group);
            operations.splice(operations.end(), group_operations);
            group.clear();
        }
        group.push_back(*it);
        last_extractor_type = extractor_type;
    }
    if (!group.empty()) {
        Independent_operations
                group_operations = operation_extractor_map_sptr->get_extractor(
                        extractor_type)->extract(reference_particle, group);
        operations.splice(operations.end(), group_operations);
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
        Operation_extractor_map_sptr const& operation_extractor_map_sptr) :
    Operator(name, "independent"), have_operations(false),
            operation_extractor_map_sptr(operation_extractor_map_sptr)
{
}

void
Independent_operator::append_slice(Lattice_element_slice_sptr const& slice_sptr)
{
    slices.push_back(slice_sptr);
}

Lattice_element_slices const&
Independent_operator::get_slices() const
{
    return slices;
}

void
Independent_operator::apply(Bunch & bunch, Operators & step_operators)
{
    if (need_update()) {
        update_operations(bunch.get_reference_particle());
    }
    for (Independent_operations::iterator it = operations.begin(); it
            != operations.end(); ++it) {
        (*it)->apply(bunch);
    }
}

void
Independent_operator::print() const
{
    Operator::print();
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        (*it)->print();
    }
}

Independent_operator::~Independent_operator()
{

}
