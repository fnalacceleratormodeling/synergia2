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
    Lattice lattice;
    Chef_lattice chef_lattice;
    Operation_extractor_map extractor_map;
    int map_order;

public:
    Lattice_simulator(Lattice & lattice, int map_order);
    void
    construct_sliced_chef_beamline(Steps const& steps);
    int
    get_map_order() const;
    void
    set_extractor(std::string const& name, Operation_extractor_sptr extractor);
    Operation_extractor_sptr
    get_extractor(std::string const& name);
    std::list<std::string >
    get_extractor_names() const;
    Operation_extractor_map &
    get_operation_extraction_map();
    double
    get_length();
    Lattice &
    get_lattice();
    Chef_lattice &
    get_chef_lattice();
    ~Lattice_simulator();
};

#endif /* LATTICE_SIMULATOR_H_ */
