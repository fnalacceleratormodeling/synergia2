#include "operation_extractor.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/lattice/chef_lattice_section.h"
#include "fast_mapping.h"
#include <cstring>

Operation_extractor::Operation_extractor(Chef_lattice_sptr chef_lattice_sptr,
        int map_order) :
    chef_lattice_sptr(chef_lattice_sptr), map_order(map_order)
{

}

Operation_extractor::Operation_extractor()
{

}

Chef_lattice_sptr &
Operation_extractor::get_chef_lattice_sptr()
{
    return chef_lattice_sptr;
}

int
Operation_extractor::get_map_order() const
{
    return map_order;
}

Operation_extractor::~Operation_extractor()
{
}

// extract_fast_mapping is a local function
Fast_mapping_operation_sptr
extract_fast_mapping(Reference_particle const& reference_particle,
        Chef_lattice_section & lattice_section, int map_order)
{
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            reference_particle, map_order);
    Particle particle = reference_particle_to_chef_particle(reference_particle);
    double mapping_length = 0.0;
    for (Chef_lattice_section::const_iterator cls_it = lattice_section.begin(); cls_it
            != lattice_section.end(); ++cls_it) {
        (*cls_it)->propagate(jet_particle);
        mapping_length += (*cls_it)->OrbitLength(particle);
        (*cls_it)->propagate(particle);
    }
    Fast_mapping fast_mapping(reference_particle, jet_particle.State(),
            mapping_length);
    Fast_mapping_operation_sptr fast_mapping_operation_sptr(
            new Fast_mapping_operation(fast_mapping));
    return fast_mapping_operation_sptr;
}

Chef_map_operation_extractor::Chef_map_operation_extractor(
        Chef_lattice_sptr chef_lattice_sptr, int map_order) :
    Operation_extractor(chef_lattice_sptr, map_order)
{
}

Chef_map_operation_extractor::Chef_map_operation_extractor()
{
}

Independent_operations
Chef_map_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Chef_lattice_section_sptr entire_section_sptr(
            new Chef_lattice_section(get_chef_lattice_sptr()));
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_lattice_section_sptr slice_section_sptr(
                get_chef_lattice_sptr()->get_chef_section_sptr(*(*les_it)));
        entire_section_sptr->extend(*slice_section_sptr);
    }
    Fast_mapping_operation_sptr fast_mapping_operation_sptr =
            extract_fast_mapping(reference_particle, *entire_section_sptr,
                    get_map_order());
    Independent_operations retval;
    retval.push_back(fast_mapping_operation_sptr);
    return retval;
}
BOOST_CLASS_EXPORT(Chef_map_operation_extractor)

Chef_propagate_operation_extractor::Chef_propagate_operation_extractor(
        Chef_lattice_sptr chef_lattice_sptr, int map_order) :
    Operation_extractor(chef_lattice_sptr, map_order)
{
}

Chef_propagate_operation_extractor::Chef_propagate_operation_extractor()
{
}

Independent_operations
Chef_propagate_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Chef_lattice_section_sptr entire_section_sptr(
            new Chef_lattice_section(get_chef_lattice_sptr()));
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_lattice_section_sptr slice_section_sptr(
                get_chef_lattice_sptr()->get_chef_section_sptr(*(*les_it)));
        entire_section_sptr->extend(*slice_section_sptr);
    }
    Chef_propagate_operation_sptr chef_propagate_operation_sptr(
            new Chef_propagate_operation(entire_section_sptr));
    Independent_operations retval;
    retval.push_back(chef_propagate_operation_sptr);
    return retval;
}
BOOST_CLASS_EXPORT(Chef_propagate_operation_extractor)

Chef_mixed_operation_extractor::Chef_mixed_operation_extractor(
        Chef_lattice_sptr chef_lattice_sptr, int map_order) :
    Operation_extractor(chef_lattice_sptr, map_order)
{
}

Chef_mixed_operation_extractor::Chef_mixed_operation_extractor()
{
}

// handle_subsection is a local function
void
handle_subsection(bool is_rf, Reference_particle const& reference_particle,
        Chef_lattice_section_sptr section_sptr,
        Independent_operations & independent_operations, int map_order)
{
    if (is_rf) {
        Chef_propagate_operation_sptr chef_propagate_operation_sptr(
                new Chef_propagate_operation(section_sptr));
        independent_operations.push_back(chef_propagate_operation_sptr);
    } else {
        Fast_mapping_operation_sptr fast_mapping_oeration_sptr =
                extract_fast_mapping(reference_particle, *section_sptr,
                        map_order);
        independent_operations.push_back(fast_mapping_oeration_sptr);
    }
}

Independent_operations
Chef_mixed_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Independent_operations retval;
    Chef_lattice_section_sptr subsection_sptr(
            new Chef_lattice_section(get_chef_lattice_sptr()));
    bool is_rf(false), last_is_rf(false);
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_lattice_section_sptr slice_section_sptr(
                get_chef_lattice_sptr()->get_chef_section_sptr(*(*les_it)));
        int index = slice_section_sptr->get_begin_index();
        for (Chef_lattice_section::const_iterator cls_it =
                slice_section_sptr->begin(); cls_it
                != slice_section_sptr->end(); ++cls_it) {
            is_rf = ((std::strcmp((*cls_it)->Type(), "rfcavity") == 0)
                    || (std::strcmp((*cls_it)->Type(), "thinrfcavity") == 0));
            if ((is_rf != last_is_rf) && (!subsection_sptr->empty())) {
                handle_subsection(last_is_rf, reference_particle,
                        subsection_sptr, retval, get_map_order());
                subsection_sptr->clear();
            }
            subsection_sptr->extend(index,index+1);
            last_is_rf = is_rf;
            ++index;
        }
    }
    if (!subsection_sptr->empty()) {
        handle_subsection(is_rf, reference_particle, subsection_sptr, retval,
                get_map_order());
    }

    return retval;
}
BOOST_CLASS_EXPORT(Chef_mixed_operation_extractor)

Operation_extractor_map::Operation_extractor_map()
{
}

void
Operation_extractor_map::set_extractor(std::string const& name,
        Operation_extractor_sptr operation_extractor)
{
    extractor_map[name] = operation_extractor;
}

Operation_extractor_sptr &
Operation_extractor_map::get_extractor(std::string const& name)
{
    return extractor_map[name];
}

std::list<std::string >
Operation_extractor_map::get_extractor_names() const
{
    std::list<std::string > retval;
    for (std::map<std::string, Operation_extractor_sptr >::const_iterator it =
            extractor_map.begin(); it != extractor_map.end(); ++it) {
        retval.push_back(it->first);
    }
    return retval;
}

Operation_extractor_map::~Operation_extractor_map()
{
}
