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
Lattice_simulator::construct_sliced_chef_beamline(Lattice_element_slices const& slices)
{
    chef_lattice_sptr->construct_sliced_beamline(slices);
}

int
Lattice_simulator::get_map_order() const
{
    return map_order;
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

