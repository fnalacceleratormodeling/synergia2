#include "communicate.h"
#include <mpi.h>

void
gather_rho(Real_scalar_field &rho)
{
    // jfa: this is the stupid way.
    void* tmp = malloc(rho.get_points().get_length() * sizeof(double));
    MPI_Allreduce(reinterpret_cast<void*>(rho.get_points().get_base_address()),
                  tmp, rho.get_points().get_length(),
                  MPI_DOUBLE, MPI_SUM, MPI_COMM_WORLD);
    memcpy(reinterpret_cast<void*>(rho.get_points().get_base_address()),
           tmp, rho.get_points().get_length()*sizeof(double));
    free(tmp);
}
