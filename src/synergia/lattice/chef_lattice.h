#ifndef CHEF_LATTICE_H_
#define CHEF_LATTICE_H_

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/chef_elements.h"

class Chef_lattice
{
    friend class Chef_lattice_tester;
private:
    Lattice_sptr lattice_sptr;
    Element_adaptor_map_sptr element_adaptor_map_sptr;
    BmlPtr beamline_sptr;
    BmlPtr sliced_beamline_sptr;
    Lattice_element_slices lattice_element_slices;
    bool have_sliced_beamline_;
    ElmPtr lattice_element_marker;
    std::map<const Lattice_element*, Chef_elements > element_map;
    std::map<const Lattice_element_slice*, Chef_elements > element_slice_map;
    double brho;

    void
    construct_beamline();
    void
    register_beamline(beamline & the_beamline);
    BmlPtr
    polish_beamline(BmlPtr beamline_sptr);
    void
    extract_element_map();
    void
    extract_element_slice_map(Lattice_element_slices const& slices);
    Chef_elements
    get_chef_elements_from_slice(Lattice_element_slice const& slice);
    void
    construct();
public:
    Chef_lattice(Lattice_sptr lattice_sptr);
    Chef_lattice(Lattice_sptr lattice_sptr,
            Element_adaptor_map_sptr element_adaptor_map_sptr);
    Element_adaptor_map_sptr
    get_element_adaptor_map_sptr();
    Chef_elements &
    get_chef_elements(Lattice_element const& lattice_element);
    Chef_elements &
    get_chef_elements(Lattice_element_slice const& lattice_element_slice);
    Lattice_element const&
    get_lattice_element(ElmPtr const& chef_element);
    Lattice_element_slice const&
    get_lattice_element_slice(ElmPtr const& chef_element);
    bool
    have_sliced_beamline() const;
    void
    construct_sliced_beamline(Lattice_element_slices const& slices);
    BmlPtr
    get_beamline_sptr();
    BmlPtr
    get_sliced_beamline_sptr();
    Lattice_element_slices const&
    get_lattice_element_slices() const;
    ~Chef_lattice();
};

typedef boost::shared_ptr<Chef_lattice > Chef_lattice_sptr;

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

#endif /* CHEF_LATTICE_H_ */
