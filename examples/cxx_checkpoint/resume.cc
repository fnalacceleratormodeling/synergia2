#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"
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
// be cleaned up properly before we call MPI_Finalize.
void
run()
{
    Propagator propagator;
    remove_serialization_directory();
    symlink_serialization_directory(Propagator::default_checkpoint_dir);
    binary_load(
            propagator,
            get_combined_path(Propagator::default_checkpoint_dir, "propagator.bina").c_str());
    unlink_serialization_directory();
    propagator.resume(Propagator::default_checkpoint_dir);
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
