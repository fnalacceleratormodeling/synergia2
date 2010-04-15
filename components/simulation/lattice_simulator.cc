#include "lattice_simulator.h"

void
Lattice_simulator::construct_extractor_map()
{
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    extractor_map.set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map.set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map.set_extractor(chef_propagate_operation_extractor_name,
            Operation_extractor_sptr(new Chef_propagate_operation_extractor(
                    chef_lattice_sptr, map_order)));
    extractor_map.set_extractor(chef_map_operation_extractor_name,
            Operation_extractor_sptr(new Chef_map_operation_extractor(
                    chef_lattice_sptr, map_order)));
}

Lattice_simulator::Lattice_simulator(Lattice_sptr const& lattice_sptr,
        int map_order) :
    lattice_sptr(lattice_sptr), chef_lattice_sptr(new Chef_lattice(
            *lattice_sptr)), map_order(map_order)
{
    construct_extractor_map();
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
    chef_lattice_sptr->construct_sliced_beamline(all_slices);
}

int
Lattice_simulator::get_map_order() const
{
    return map_order;
}

void
Lattice_simulator::set_extractor(std::string const& name,
        Operation_extractor_sptr extractor)
{
    extractor_map.set_extractor(name, extractor);
}

Operation_extractor_sptr
Lattice_simulator::get_extractor(std::string const& name)
{
    return extractor_map.get_extractor(name);
}

std::list<std::string >
Lattice_simulator::get_extractor_names() const
{
}

Operation_extractor_map &
Lattice_simulator::get_operation_extraction_map()
{
    return extractor_map;
}

double
Lattice_simulator::get_length()
{
}

Lattice_sptr &
Lattice_simulator::get_lattice_sptr()
{
    return lattice_sptr;
}

Chef_lattice_sptr &
Lattice_simulator::get_chef_lattice_sptr()
{
    return chef_lattice_sptr;
}

Lattice_simulator::~Lattice_simulator()
{
}

