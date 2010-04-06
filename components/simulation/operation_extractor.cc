#include "operation_extractor.h"
#include "components/lattice/chef_utils.h"
#include "fast_mapping.h"

// extract_fast_mapping is a local function
Fast_mapping_operation_sptr
extract_fast_mapping(Reference_particle const& reference_particle,
        Chef_elements & chef_elements)
{
    JetParticle jet_particle = reference_particle_to_chef_jet_particle(
            reference_particle);
    for (Chef_elements::const_iterator ce_it = chef_elements.begin(); ce_it
            != chef_elements.end(); ++ce_it) {
        (*ce_it)->propagate(jet_particle);
    }
    Fast_mapping fast_mapping(reference_particle, jet_particle.State());
    Fast_mapping_operation_sptr fast_mapping_operation_sptr(
            new Fast_mapping_operation(fast_mapping));
    return fast_mapping_operation_sptr;
}

Independent_operations
Chef_map_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices, Chef_lattice & chef_lattice)
{
    Chef_elements all_chef_elements;
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_elements slice_elements = chef_lattice.get_chef_elements(
                *(*les_it));
        for (Chef_elements::const_iterator ce_it = slice_elements.begin(); ce_it
                != slice_elements.end(); ++ce_it) {
            all_chef_elements.push_back(*ce_it);
        }
    }
    Fast_mapping_operation_sptr fast_mapping_operation_sptr =
            extract_fast_mapping(reference_particle, all_chef_elements);
    Independent_operations retval;
    retval.push_back(fast_mapping_operation_sptr);
    return retval;
}

Independent_operations
Chef_propagate_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices, Chef_lattice & chef_lattice)
{
    Chef_elements chef_elements;
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_elements slice_elements = chef_lattice.get_chef_elements(
                *(*les_it));
        for (Chef_elements::const_iterator ce_it = slice_elements.begin(); ce_it
                != slice_elements.end(); ++ce_it) {
            chef_elements.push_back(*ce_it);
        }
    }
    Chef_propagate_operation_sptr chef_propagate_operation_sptr(
            new Chef_propagate_operation(chef_elements));
    Independent_operations retval;
    retval.push_back(chef_propagate_operation_sptr);
    return retval;
}

// handle_group is a local function
void
handle_group(bool is_rf, Reference_particle const& reference_particle,
        Chef_elements & group, Independent_operations & independent_operations)
{
    if (is_rf) {
        Chef_propagate_operation_sptr chef_propagate_operation_sptr(
                new Chef_propagate_operation(group));
        independent_operations.push_back(chef_propagate_operation_sptr);
    } else {
        Fast_mapping_operation_sptr fast_mapping_oeration_sptr =
                extract_fast_mapping(reference_particle, group);
        independent_operations.push_back(fast_mapping_oeration_sptr);
    }
}

Independent_operations
Mixed_chef_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices, Chef_lattice & chef_lattice)
{
    Independent_operations retval;
    Chef_elements group;
    bool is_rf(false), last_is_rf(false);
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_elements slice_elements = chef_lattice.get_chef_elements(
                *(*les_it));
        for (Chef_elements::const_iterator ce_it = slice_elements.begin(); ce_it
                != slice_elements.end(); ++ce_it) {
            is_rf = (((*ce_it)->Type() == "rfcavity") || ((*ce_it)->Type()
                    == "thinrfcavity"));
            if ((is_rf != last_is_rf) && (!group.empty())) {
                handle_group(last_is_rf, reference_particle, group, retval);
                group.clear();
            }
            group.push_back(*ce_it);
            last_is_rf = is_rf;
        }
    }
    if (!group.empty()) {
        handle_group(is_rf, reference_particle, group, retval);
    }

    return retval;
}
