#ifndef CHEF_LATTICE_H_
#define CHEF_LATTICE_H_

#include "components/lattice/lattice.h"
#include <beamline/beamline.h>

typedef std::list<ElmPtr > Chef_elements;
typedef Chef_elements
(*lattice_element_to_chef_fn)(Lattice_element const&, double);
typedef std::map<std::string, lattice_element_to_chef_fn >
        Lattice_element_to_chef_fn_map;

class Chef_lattice
{
private:
    Lattice *lattice_ptr;
    BmlPtr beamline_ptr;
    double brho;

    beamline
    construct_raw_lattice(Lattice_element_to_chef_fn_map const& map);
    void
    polish_lattice(beamline const& raw_beamline);
    void
    construct(Lattice_element_to_chef_fn_map const& map);
public:
    Chef_lattice(Lattice & lattice);
    Chef_lattice(Lattice & lattice, Lattice_element_to_chef_fn_map const& map);
    BmlPtr
    get_beamline_ptr();
    ~Chef_lattice();
};

Particle
reference_particle_to_chef_particle(
        Reference_particle const& reference_particle);

Lattice_element_to_chef_fn_map
get_standard_lattice_element_to_chef_fn_map();

Chef_elements
lattice_element_to_chef_marker(Lattice_element const& lattice_element,
        double brho);
Chef_elements
lattice_element_to_chef_drift(Lattice_element const& lattice_element,
        double brho);
Chef_elements
lattice_element_to_chef_quadrupole(Lattice_element const& lattice_element,
        double brho);

#endif /* CHEF_LATTICE_H_ */
