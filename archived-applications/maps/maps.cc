#include <iostream>
#include <mpi.h>
#include <stdexcept>

#include "synergia/foundation/trigon_particle.h"
#include "synergia/lattice/chef_utils.h"
#include "synergia/lattice/lattice.h"
#include "synergia/lattice/lattice_element_slice.h"
#include "synergia/lattice/madx_reader.h"
#include "synergia/libFF/ff_element_map.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    MadX_reader reader;
    Lattice_sptr lattice_sptr = reader.get_lattice_sptr("fodo", "fodo.madx");
    //    for(auto it = lattice_sptr->get_elements().begin();
    //        it != lattice_sptr->get_elements().end(); ++it){
    //        std::cout << (*it)->get_name() << std::endl;
    //    }

    construct_big_giant_global_ff_element_map();
    FF_element_map& element_map(the_big_giant_global_ff_element_map);

    JetParticle jet_particle(reference_particle_to_chef_jet_particle(
        lattice_sptr->get_reference_particle(),
        Trigon_particle_t::Component_t::power()));
    for (auto&& element_sptr : lattice_sptr->get_elements()) {
        Lattice_element_slice slice(element_sptr);
        std::cout << element_sptr->get_name() << ": "
                  << element_sptr->get_type() << std::endl;
        element_map.get_element_type(element_sptr->get_type())
            ->apply(element_sptr, jet_particle);
        std::cout << "jet_particle.jacobian() = \n";
        std::cout << jet_particle.State().Jacobian();
    }
}

int
main(int argc, char** argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
