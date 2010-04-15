#include "operation_extractor.h"
#include "components/lattice/chef_utils.h"
#include "fast_mapping.h"
#include <cstring>

Operation_extractor::Operation_extractor(
        Chef_lattice_sptr const& chef_lattice_sptr, int map_order) :
    chef_lattice_sptr(chef_lattice_sptr), map_order(map_order)
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
        Chef_elements & chef_elements, int map_order)
{
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            reference_particle, map_order);
    for (Chef_elements::const_iterator ce_it = chef_elements.begin(); ce_it
            != chef_elements.end(); ++ce_it) {
        (*ce_it)->propagate(jet_particle);
    }
    Fast_mapping fast_mapping(reference_particle, jet_particle.State());
    Fast_mapping_operation_sptr fast_mapping_operation_sptr(
            new Fast_mapping_operation(fast_mapping));
    return fast_mapping_operation_sptr;
}

Chef_map_operation_extractor::Chef_map_operation_extractor(
        Chef_lattice_sptr const& chef_lattice_sptr, int map_order) :
    Operation_extractor(chef_lattice_sptr, map_order)
{
}

Independent_operations
Chef_map_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Chef_elements all_chef_elements;
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_elements slice_elements =
                get_chef_lattice_sptr()->get_chef_elements(*(*les_it));
        all_chef_elements.splice(all_chef_elements.end(), slice_elements);
    }
    Fast_mapping_operation_sptr fast_mapping_operation_sptr =
            extract_fast_mapping(reference_particle, all_chef_elements,
                    get_map_order());
    Independent_operations retval;
    retval.push_back(fast_mapping_operation_sptr);
    return retval;
}

Chef_propagate_operation_extractor::Chef_propagate_operation_extractor(
        Chef_lattice_sptr const& chef_lattice_sptr, int map_order) :
    Operation_extractor(chef_lattice_sptr, map_order)
{
}

Independent_operations
Chef_propagate_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Chef_elements chef_elements;
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_elements slice_elements =
                get_chef_lattice_sptr()->get_chef_elements(*(*les_it));
        chef_elements.splice(chef_elements.end(), slice_elements);
    }
    Chef_propagate_operation_sptr chef_propagate_operation_sptr(
            new Chef_propagate_operation(chef_elements));
    Independent_operations retval;
    retval.push_back(chef_propagate_operation_sptr);
    return retval;
}

Chef_mixed_operation_extractor::Chef_mixed_operation_extractor(
        Chef_lattice_sptr const& chef_lattice_sptr, int map_order) :
    Operation_extractor(chef_lattice_sptr, map_order)
{
}

// handle_group is a local function
void
handle_group(bool is_rf, Reference_particle const& reference_particle,
        Chef_elements & group, Independent_operations & independent_operations,
        int map_order)
{
    if (is_rf) {
        Chef_propagate_operation_sptr chef_propagate_operation_sptr(
                new Chef_propagate_operation(group));
        independent_operations.push_back(chef_propagate_operation_sptr);
    } else {
        Fast_mapping_operation_sptr fast_mapping_oeration_sptr =
                extract_fast_mapping(reference_particle, group, map_order);
        independent_operations.push_back(fast_mapping_oeration_sptr);
    }
}

Independent_operations
Chef_mixed_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Independent_operations retval;
    Chef_elements group;
    bool is_rf(false), last_is_rf(false);
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_elements slice_elements =
                get_chef_lattice_sptr()->get_chef_elements(*(*les_it));
        for (Chef_elements::const_iterator ce_it = slice_elements.begin(); ce_it
                != slice_elements.end(); ++ce_it) {
            is_rf = ((std::strcmp((*ce_it)->Type(), "rfcavity") == 0)
                    || (std::strcmp((*ce_it)->Type(), "thinrfcavity") == 0));
            if ((is_rf != last_is_rf) && (!group.empty())) {
                handle_group(last_is_rf, reference_particle, group, retval,
                        get_map_order());
                group.clear();
            }
            group.push_back(*ce_it);
            last_is_rf = is_rf;
        }
    }
    if (!group.empty()) {
        handle_group(is_rf, reference_particle, group, retval, get_map_order());
    }

    return retval;
}

Operation_extractor_map::Operation_extractor_map()
{
}

void
Operation_extractor_map::set_extractor(std::string const& name,
        Operation_extractor_sptr const& operation_extractor)
{
    extractor_map[name] = operation_extractor;
}

Operation_extractor_sptr &
Operation_extractor_map::get_extractor(std::string const& name)
{
    return extractor_map[name];
}

Operation_extractor_map::~Operation_extractor_map()
{
}
