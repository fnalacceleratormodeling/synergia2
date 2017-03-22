#include <iostream>
#include <mpi.h>
#include <stdexcept>

#include "synergia/foundation/trigon_particle.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/libFF/ff_element_map.h"
#include "synergia/utils/multi_array_print.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    MadX_reader reader;
    Lattice_sptr lattice_sptr = reader.get_lattice_sptr("model", "foborodobo128.madx");

    construct_big_giant_global_ff_element_map();
    FF_element_map& element_map(the_big_giant_global_ff_element_map);

    JetParticle jet_particle(reference_particle_to_chef_jet_particle(
        lattice_sptr->get_reference_particle(),
        Trigon_particle_t::Component_t::power()));

    double t0 = MPI_Wtime();
    const int num_reps = 5;
    for (int i = 0; i < num_reps; ++i) {
        for (auto&& element_sptr : lattice_sptr->get_elements()) {
//            std::cout << element_sptr->get_name() << ": " <<
//                         element_sptr->get_type() << std::endl;
            Lattice_element_slice slice(element_sptr);
            element_map.get_element_type(element_sptr->get_type())
                ->apply(element_sptr, jet_particle);
        }
    }
    double t1 = MPI_Wtime();
    std::cout << "chef time: " << t1 - t0 << std::endl;
    std::cout << "jet_particle.jacobian() = \n";
    std::cout << jet_particle.State().Jacobian();
//    jet_particle.State().printCoeffs();

    Trigon_particle_t trigon_particle(lattice_sptr->get_reference_particle());

    t0 = MPI_Wtime();
    for (int i = 0; i < num_reps; ++i) {
        for (auto&& element_sptr : lattice_sptr->get_elements()) {
            Lattice_element_slice slice(element_sptr);
            element_map.get_element_type(element_sptr->get_type())
                ->apply(element_sptr, trigon_particle);
        }
    }
    t1 = MPI_Wtime();
    std::cout << "trigon time: " << t1 - t0 << std::endl;
    std::cout << "trigon_particle.jacobian() = \n";
    multi_array_print(trigon_particle.get_jacobian(), "jacobian");
}

int
main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
