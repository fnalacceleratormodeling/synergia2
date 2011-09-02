#include "lattice_simulator.h"

void
Lattice_simulator::construct_extractor_map()
{
    Operation_extractor_sptr chef_mixed_operation_extractor(
            new Chef_mixed_operation_extractor(chef_lattice_sptr, map_order));

    extractor_map_sptr->set_extractor(default_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(chef_mixed_operation_extractor_name,
            chef_mixed_operation_extractor);
    extractor_map_sptr->set_extractor(
            chef_propagate_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_propagate_operation_extractor(chef_lattice_sptr,
                            map_order)));
    extractor_map_sptr->set_extractor(
            chef_map_operation_extractor_name,
            Operation_extractor_sptr(
                    new Chef_map_operation_extractor(chef_lattice_sptr,
                            map_order)));
}

Lattice_simulator::Lattice_simulator(Lattice_sptr lattice_sptr, int map_order) :
    lattice_sptr(lattice_sptr),
            chef_lattice_sptr(new Chef_lattice(lattice_sptr)),
            extractor_map_sptr(new Operation_extractor_map),
            map_order(map_order), have_slices(false)
{
    construct_extractor_map();
}

void
Lattice_simulator::set_slices(Lattice_element_slices const& slices)
{
    this->slices = slices;
    have_slices = true;
    construct_sliced_chef_beamline();
}

void
Lattice_simulator::construct_sliced_chef_beamline()
{
    if (!have_slices) {
        throw std::runtime_error(
                "Lattice_simulator::construct_sliced_chef_beamline called before set_slices");
    }
    chef_lattice_sptr->construct_sliced_beamline(slices);
}

int
Lattice_simulator::get_map_order() const
{
    return map_order;
}

Operation_extractor_map_sptr
Lattice_simulator::get_operation_extractor_map_sptr()
{
    return extractor_map_sptr;
}

Lattice_sptr
Lattice_simulator::get_lattice_sptr()
{
    return lattice_sptr;
}

Chef_lattice_sptr
Lattice_simulator::get_chef_lattice_sptr()
{
    return chef_lattice_sptr;
}

void
Lattice_simulator::update()
{
    chef_lattice_sptr = Chef_lattice_sptr(new Chef_lattice(lattice_sptr));
    construct_extractor_map();
    if (have_slices) {
        construct_sliced_chef_beamline();
    }
}

Lattice_simulator::~Lattice_simulator()
{
}

