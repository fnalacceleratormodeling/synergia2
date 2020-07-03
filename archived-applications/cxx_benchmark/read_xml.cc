#include <iostream>
#include <stdexcept>

#include "synergia/lattice/lattice.h"
#include "synergia/utils/serialization.h"

#include "synergia/utils/synergia_omp.h"

// We put the actual code in a separate function so that shared_ptr's can
// be cleanup up properly before we call MPI_Finalize.
void
run()
{
    MArray2d covariances;
    xml_load(covariances, "covariance.xml");
    std::cout.precision(16);
    for (int i=0; i<6; ++i) {
        if (i != 0) {
            std::cout << std::endl;
        }
        for (int j=0; j<6; ++j) {
            if (j != 0) {
                std::cout << " ";
            }
            std::cout << covariances[i][j];
        }
    }
    std::cout << std::endl;
}

int
main(int argc, char **argv)
{
    MPI_Init(&argc, &argv);
    run();
    MPI_Finalize();
    return 0;
}
