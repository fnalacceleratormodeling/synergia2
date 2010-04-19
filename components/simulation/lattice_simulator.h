#ifndef LATTICE_SIMULATOR_H_
#define LATTICE_SIMULATOR_H_

#include "components/lattice/lattice.h"
#include "components/lattice/chef_lattice.h"
#include "components/simulation/operation_extractor.h"
#include "components/simulation/step.h"
#include <string>

class Lattice_simulator
{
private:
    Lattice_sptr lattice_sptr;
    Chef_lattice_sptr chef_lattice_sptr;
    Operation_extractor_map_sptr extractor_map_sptr;
    int map_order;

    void
    construct_extractor_map();
public:
    Lattice_simulator(Lattice_sptr const& lattice, int map_order);
    void
    construct_sliced_chef_beamline(Lattice_element_slices const& slices);
    int
    get_map_order() const;
    Operation_extractor_map_sptr &
    get_operation_extractor_map_sptr();
    Lattice_sptr &
    get_lattice_sptr();
    Chef_lattice_sptr &
    get_chef_lattice_sptr();
    ~Lattice_simulator();
};

#endif /* LATTICE_SIMULATOR_H_ */
