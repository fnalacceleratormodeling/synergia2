#include <iostream>
#include "synergia/utils/simple_timer.h"
#include "operator.h"
#include "aperture.h"

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

Collective_operator::~Collective_operator()
{
}

Dummy_collective_operator::Dummy_collective_operator(std::string const& name) :
    Collective_operator(name)
{
}

void
Dummy_collective_operator::apply(Bunch & bunch, double time_step, Step & step)
{
}

Dummy_collective_operator::~Dummy_collective_operator()
{
}

void
Independent_operator::update_operations(
        Reference_particle const& reference_particle)
{
    operations.clear();
    operations_revisions.clear();

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
        operations_revisions.push_back(
                (*it)->get_lattice_element().get_revision());
    }
    if (!group.empty()) {
        Independent_operations
                group_operations = operation_extractor_map_sptr->get_extractor(
                        extractor_type)->extract(reference_particle, group);
        operations.splice(operations.end(), group_operations);
    }
    have_operations = true;
    operations_reference_particle = reference_particle;
}

bool
Independent_operator::need_update(Reference_particle const& reference_particle)
{
    const double reference_particle_tolerance = 1.0e-8;
    bool retval;
    if (have_operations) {
        retval = false;
        if (reference_particle.equal(operations_reference_particle,
                reference_particle_tolerance)) {
            std::list<long int >::const_iterator rev_it =
                    operations_revisions.begin();
            for (Lattice_element_slices::const_iterator it = slices.begin(); it
                    != slices.end(); ++it) {
                long int cached_revision = (*rev_it);
                long int revision = (*it)->get_lattice_element().get_revision();
                if (revision != cached_revision) {
                    retval = true;
                }
            }
        } else {
            retval = true;
        }
    } else {
        retval = true;
    }
    return retval;
}

Independent_operator::Independent_operator(std::string const& name,
        Operation_extractor_map_sptr operation_extractor_map_sptr) :
    Operator(name, "independent"),
            operation_extractor_map_sptr(operation_extractor_map_sptr),
            have_operations(false)
{
}

void
Independent_operator::append_slice(Lattice_element_slice_sptr slice_sptr)
{
    slices.push_back(slice_sptr);
    have_operations = false;
}

Lattice_element_slices const&
Independent_operator::get_slices() const
{
    return slices;
}

void
Independent_operator::apply(Bunch & bunch, double time_step, Step & step)
{
    if (need_update(bunch.get_reference_particle())) {
        update_operations(bunch.get_reference_particle());
    }
    double t;
    t = simple_timer_current();
    for (Independent_operations::iterator it = operations.begin(); it
            != operations.end(); ++it) {
        (*it)->apply(bunch);
    }
    t = simple_timer_show(t, "independent-operation-apply");
    apply_circular_aperture(bunch, slices);
    t = simple_timer_show(t, "independent-operation-aperture");
}

void
Independent_operator::print() const
{
    Operator::print();
    int count = 0;
    for (Lattice_element_slices::const_iterator it = slices.begin(); it
            != slices.end(); ++it) {
        ++count;
        std::cout << "slice " << count << ": ";
        (*it)->print();
    }
}

Independent_operator::~Independent_operator()
{

}
