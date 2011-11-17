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
    double arc_length;
};

class Lattice_simulator
{
private:
    Lattice_sptr lattice_sptr;
    Lattice_element_slices slices;
    bool have_slices;
    Chef_lattice_sptr chef_lattice_sptr;
    Operation_extractor_map_sptr extractor_map_sptr;
    int map_order;
    double bucket_length;
    void
    construct_extractor_map();
    void
    construct_sliced_chef_beamline();
    bool have_element_lattice_functions;
    bool have_slice_lattice_functions;
    double horizontal_tune, vertical_tune;
    bool have_tunes;
    std::map<Lattice_element const*, Lattice_functions >
            lattice_functions_element_map;
    std::map<Lattice_element_slice const*, Lattice_functions >
            lattice_functions_slice_map;
    void
    get_tunes();
public:
    Lattice_simulator(Lattice_sptr lattice, int map_order);
    void
    set_slices(Lattice_element_slices const& slices);
    int
    get_map_order() const;
    void
    set_bucket_length();
    /// bucket length is in z_lab frame
    double
    get_bucket_length();
    int
    get_number_buckets();
    Operation_extractor_map_sptr
    get_operation_extractor_map_sptr();
    Lattice_sptr
    get_lattice_sptr();
    Chef_lattice_sptr
    get_chef_lattice_sptr();
    void
    update();
    void
    calculate_element_lattice_functions();
    void
    calculate_slice_lattice_functions();
    Lattice_functions const&
    get_lattice_functions(Lattice_element & lattice_element);
    Lattice_functions const&
    get_lattice_functions(Lattice_element_slice & lattice_element_slice);
    double
    get_horizontal_tune();
    double
    get_vertical_tune();
    void
    adjust_tunes(double horizontal_tune, double vertical_tune,
            Lattice_elements const& horizontal_correctors,
            Lattice_elements const& vertical_correctors,
            double tolerance = 1.0e-6);
    ~Lattice_simulator();
};

#endif /* LATTICE_SIMULATOR_H_ */
