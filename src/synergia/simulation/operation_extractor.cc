
#include "synergia/simulation/operation_extractor.h"

#include "synergia/simulation/independent_operation.h"
#include "synergia/simulation/aperture_operation.h"

namespace
{

#if 0
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
    Fast_mapping fast_mapping(jet_particle.State(),
            mapping_length);
    Fast_mapping_operation_sptr fast_mapping_operation_sptr(
            new Fast_mapping_operation(fast_mapping));
    return fast_mapping_operation_sptr;
}

Independent_operations
Chef_map_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Chef_lattice_section entire_section(get_chef_lattice_sptr());
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_lattice_section slice_section(
                *(get_chef_lattice_sptr()->get_chef_section_sptr(get_chef_lattice_sptr(),
                        *(*les_it))));
        entire_section.extend(slice_section);
    }
    Fast_mapping_operation_sptr fast_mapping_operation_sptr =
            extract_fast_mapping(reference_particle, entire_section,
                    get_map_order());
    Independent_operations retval;
    retval.push_back(fast_mapping_operation_sptr);
    return retval;
}



Independent_operations
Chef_propagate_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Chef_lattice_section entire_section(get_chef_lattice_sptr());
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_lattice_section slice_section(
                *(get_chef_lattice_sptr()->get_chef_section_sptr(get_chef_lattice_sptr(),
                        *(*les_it))));
        entire_section.extend(slice_section);
    }
    Chef_propagate_operation_sptr chef_propagate_operation_sptr(
            new Chef_propagate_operation(entire_section));
    Independent_operations retval;
    retval.push_back(chef_propagate_operation_sptr);
    return retval;
}


// handle_subsection is a local function
void
handle_subsection(bool is_rf, Reference_particle const& reference_particle,
        Chef_lattice_section & section,
        Independent_operations & independent_operations, int map_order)
{
    if (is_rf) {
        Chef_propagate_operation_sptr chef_propagate_operation_sptr(
                new Chef_propagate_operation(section));
        independent_operations.push_back(chef_propagate_operation_sptr);
    } else {
        Fast_mapping_operation_sptr fast_mapping_operation_sptr =
                extract_fast_mapping(reference_particle, section,
                        map_order);
        independent_operations.push_back(fast_mapping_operation_sptr);
    }
}

Independent_operations
Chef_mixed_operation_extractor::extract(
        Reference_particle const& reference_particle,
        Lattice_element_slices const& slices)
{
    Independent_operations retval;
    Chef_lattice_section subsection(get_chef_lattice_sptr());
    bool is_rf(false), last_is_rf(false);
    for (Lattice_element_slices::const_iterator les_it = slices.begin(); les_it
            != slices.end(); ++les_it) {
        Chef_lattice_section slice_section(
                *(get_chef_lattice_sptr()->get_chef_section_sptr(get_chef_lattice_sptr(),
                        *(*les_it))));
        int index = slice_section.get_begin_index();
        for (Chef_lattice_section::const_iterator cls_it =
                slice_section.begin(); cls_it
                != slice_section.end(); ++cls_it) {
            is_rf = ((std::strcmp((*cls_it)->Type(), "rfcavity") == 0)
                    || (std::strcmp((*cls_it)->Type(), "thinrfcavity") == 0));
            if ((is_rf != last_is_rf) && (!subsection.empty())) {
                handle_subsection(last_is_rf, reference_particle,
                        subsection, retval, get_map_order());
                subsection.clear();
            }
            subsection.extend(index,index+1);
            last_is_rf = is_rf;
            ++index;
        }
    }
    if (!subsection.empty()) {
        handle_subsection(is_rf, reference_particle, subsection, retval,
                get_map_order());
    }

    return retval;
}

#endif

    void
    chef_map_operation_extract(
            Lattice const & lattice,
            std::vector<Lattice_element_slice> const & slices,
            std::vector<std::unique_ptr<Independent_operation>> & operations )
    {
    }

    void
    chef_propagator_operation_extract(
            Lattice const & lattice,
            std::vector<Lattice_element_slice> const & slices,
            std::vector<std::unique_ptr<Independent_operation>> & operations )
    {
    }


    void
    libff_operation_extract(
            Lattice const & lattice,
            std::vector<Lattice_element_slice> const & slices,
            std::vector<std::unique_ptr<Independent_operation>> & operations )
    {
        operations.emplace_back(std::make_unique<LibFF_operation>(slices));
    }

} // namespace





void
extract_independent_operations(
        std::string const & extractor_type,
        Lattice const & lattice,
        std::vector<Lattice_element_slice> const & slices,
        std::vector<std::unique_ptr<Independent_operation>> & operations )
{
    if (extractor_type == "chef_map")
    {
        chef_map_operation_extract(lattice, slices, operations);
    }
    else if (extractor_type == "chef_propagator")
    {
        chef_propagator_operation_extract(lattice, slices, operations);
    }
#if 0
    else if (extractor_type == "chef_mixed")
    {
        return chef_mixed_operation_extract(lattice, slices);
    }
#endif
    else if (extractor_type == "libff" || extractor_type == "default")
    {
        libff_operation_extract(lattice, slices, operations);
    }
    else
    {
        throw std::runtime_error("unknown extractor_type");
    }
}


std::unique_ptr<Independent_operation>
extract_aperture_operation(
        std::string const & aperture_type,
        Lattice_element_slice const & slice)
{
    if (aperture_type == "finite_aperture")
    {
        return std::make_unique<Aperture_operation<Finite_aperture>>(slice);
    }
    else
    {
        throw std::runtime_error("unknown aperture_type");
    }
}




