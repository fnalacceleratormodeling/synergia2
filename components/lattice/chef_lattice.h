#ifndef CHEF_LATTICE_H_
#define CHEF_LATTICE_H_

#include "components/lattice/lattice.h"
#include "components/lattice/lattice_element_slice.h"
#include "components/lattice/chef_elements.h"

class Chef_lattice
{
private:
    Lattice *lattice_ptr;
    BmlPtr beamline_sptr;
    BmlPtr sliced_beamline_sptr;
    bool have_sliced_beamline;
    ElmPtr lattice_element_marker;
    std::map<const Lattice_element*, Chef_elements > element_map;
    std::map<const Lattice_element_slice*, Chef_elements > element_slice_map;
    double brho;

    beamline
    construct_raw_beamline();
    void
    register_beamline(beamline & the_beamline);
    void
    polish_raw_beamline(beamline const& raw_beamlinee);
    void
    extract_element_map();
    Chef_elements
    get_chef_elements_from_slice(Lattice_element_slice const& slice);
public:
    Chef_lattice(Lattice & lattice);
    Chef_elements &
    get_chef_elements(Lattice_element const& lattice_element);
    Chef_elements &
    get_chef_elements(Lattice_element_slice const& lattice_element_slice);
    void
    construct_sliced_beamline(Lattice_element_slices const& slices);
    BmlPtr
    get_beamline_sptr();
    BmlPtr
    get_sliced_beamline_sptr();
    ~Chef_lattice();
};

typedef boost::shared_ptr<Chef_lattice > Chef_lattice_sptr;

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

#endif /* CHEF_LATTICE_H_ */
