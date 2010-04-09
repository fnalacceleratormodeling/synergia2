#include "lattice_simulator.h"

Lattice_simulator::Lattice_simulator(Lattice & lattice, int map_order) :
    lattice(lattice), chef_lattice(lattice)
{
}

void
Lattice_simulator::construct_sliced_chef_beamline(Steps const& steps)
{
    Lattice_element_slices all_slices;
    for (Steps::const_iterator s_it = steps.begin(); s_it != steps.end(); ++s_it) {
        for (Operators::const_iterator o_it = (*s_it)->get_operators().begin(); o_it
                != (*s_it)->get_operators().end(); ++o_it) {
            if ((*o_it)->get_type() == "independent") {
                Lattice_element_slices
                        element_slices(boost::static_pointer_cast<
                                Independent_operator >(*o_it)->get_slices());
                all_slices.splice(all_slices.end(), element_slices);
            }
        }
    }
    chef_lattice.construct_sliced_beamline(all_slices);
}

int
Lattice_simulator::get_map_order() const
{
}

void
Lattice_simulator::set_extractor(std::string const& name,
        Operation_extractor_sptr extractor)
{
}

Operation_extractor_sptr
Lattice_simulator::get_extractor(std::string const& name)
{
}

std::list<std::string >
Lattice_simulator::get_extractor_names() const
{
}

Operation_extractor_map &
Lattice_simulator::get_operation_extraction_map()
{
}

double
Lattice_simulator::get_length()
{
}

Lattice &
Lattice_simulator::get_lattice()
{
    return lattice;
}

Chef_lattice &
Lattice_simulator::get_chef_lattice()
{
    return chef_lattice;
}

Lattice_simulator::~Lattice_simulator()
{
}

