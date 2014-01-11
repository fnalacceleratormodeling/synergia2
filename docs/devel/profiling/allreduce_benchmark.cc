#include <iostream>
#include <mpi.h>
#include "synergia/utils/commxx.h"
#include "synergia/utils/multi_array_typedefs.h"

#ifndef USE_SIMPLE_TIMER
#define USE_SIMPLE_TIMER true
#endif // USE_SIMPLE_TIMER

#include "synergia/utils/simple_timer.h"

#include "allreduce_benchmark_options.h"

int
main(int argc, char **argv)
{

    MPI_Init(&argc, &argv);
    double t = simple_timer_current();
    Allreduce_benchmark_options opts(argc, argv);

    Commxx comm;
    std::cout << "hello world from rank " << comm.get_rank() << "/"
            << comm.get_size() << "\n";

    MArray3d array(boost::extents[opts.gridnum][opts.gridnum][opts.gridnum]);
    int rank = comm.get_rank();
    for (int i = 0; i < opts.gridnum; ++i) {
        for (int j = 0; j < opts.gridnum; ++j) {
            for (int k = 0; k < opts.gridnum; ++k) {
                array[i][j][k] = rank * 1e9 + 1e6 * i + 1e3 * j + k;
            }
        }
    }

    t = simple_timer_show(t, "startup");
    for (int iter = 0; iter < opts.iterations; ++iter) {
        MPI_Allreduce(MPI_IN_PLACE, (void*) array.origin(),
                array.num_elements(), MPI_DOUBLE, MPI_SUM, comm.get());
        t = simple_timer_show(t, "allreduce");
    }

    MPI_Finalize();
    return 0;
}
