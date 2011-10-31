#ifndef LATTICE_SIMULATOR_H_
#define LATTICE_SIMULATOR_H_

#include "synergia/lattice/lattice.h"
#include "synergia/lattice/chef_lattice.h"
#include "synergia/simulation/operation_extractor.h"
#include "synergia/simulation/step.h"
#include <physics_toolkit/LattFuncSage.h>
#include <string>

struct Lattice_functions
{
    Lattice_functions();
    Lattice_functions(LattFuncSage::lattFunc const& latt_func);
    double alpha_x, alpha_y;
    double beta_x, beta_y;
    double psi_x, psi_y;
    double D_x, D_y;
    double Dprime_x, Dprime_y;
};

class Lattice_simulator
{
private:
    Lattice_sptr lattice_sptr;
    Chef_lattice_sptr chef_lattice_sptr;
    Operation_extractor_map_sptr extractor_map_sptr;
    int map_order;
    bool have_element_lattice_functions;
    bool have_slice_lattice_functions;
    std::map<Lattice_element*, Lattice_functions >
            lattice_functions_element_map;
    std::map<Lattice_element_slice*, Lattice_functions >
            lattice_functions_slice_map;
    void
    construct_extractor_map();
public:
    Lattice_simulator(Lattice_sptr lattice, int map_order);
    void
    construct_sliced_chef_beamline(Lattice_element_slices const& slices);
    int
    get_map_order() const;
    Operation_extractor_map_sptr
    get_operation_extractor_map_sptr();
    Lattice_sptr
    get_lattice_sptr();
    Chef_lattice_sptr
    get_chef_lattice_sptr();
    void
    calculate_element_lattice_functions();
    void
    calculate_slice_lattice_functions();
    Lattice_functions const&
    get_lattice_functions(Lattice_element & lattice_element);
    Lattice_functions const&
    get_lattice_functions(Lattice_element_slice & lattice_element_slice);
    ~Lattice_simulator();
};

#endif /* LATTICE_SIMULATOR_H_ */
