#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization_files.h"
#include "synergia/simulation/operator.h"
#include "synergia/simulation/lattice_simulator.h"
#include "synergia/simulation/stepper.h"
#include "synergia/simulation/propagator.h"
#include "synergia/bunch/bunch.h"
#include "synergia/foundation/distribution.h"
#include "synergia/bunch/populate.h"
#include "synergia/bunch/diagnostics.h"
#include "synergia/collective/space_charge_3d_open_hockney.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    std::vector<int > grid_shape(3);
    grid_shape[0] = 32;
    grid_shape[1] = 32;
    grid_shape[2] = 256;
    const int part_per_cell = 10;
    const int num_macro_particles = grid_shape[0] * grid_shape[1]
            * grid_shape[2] * part_per_cell;
    const int seed = 4;
    const double num_real_particles = 1e13;
    const int num_steps = 8;
    const int num_turns = 4;
    const int map_order = 3;
    const double trans_emit = 1e-6;

    Lattice_sptr lattice_sptr(new Lattice());
    try {
        xml_load(*lattice_sptr, "cxx_lattice.xml");
    }
    catch (std::runtime_error) {
        std::cerr << "normal_form_example: failed to find cxx_lattice.xml\n";
        std::cerr << "Run normal_form_example.py to generate cxx_lattice.xml\n";
        exit(1);
    }


    Lattice_simulator lattice_simulator(lattice_sptr, map_order);

    MArray2d map = lattice_simulator.get_linear_one_turn_map();

    std::cout << "First try One turn map: " << std::endl;
    std::cout << std::setprecision(15) << std::endl;
    for (int i=0; i<6; ++i) {
      for (int j=0; j<6; ++j) {
	std::cout << map[i][j] << "    ";
      }
      std::cout << endl;
    }
    
    MArray2d map2 = lattice_simulator.get_linear_one_turn_map();

    std::cout << "Second try One turn map: " << std::endl;
    std::cout << std::setprecision(10) << std::endl;
    for (int i=0; i<6; ++i) {
      for (int j=0; j<6; ++j) {
	std::cout << map2[i][j] << "    ";
      }
      std::cout << endl;
    }
}

int
main(int argc, char **argv)
{
    run();
    return 0;
}
