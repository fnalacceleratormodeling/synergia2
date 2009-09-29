#ifndef PARALLEL_UTILS_H_
#define PARALLEL_UTILS_H_

#include "mpi.h"
#include <vector>
#include <stdexcept>

int
mpi_get_rank(const MPI_Comm &comm);

int
mpi_get_size(MPI_Comm &comm);

void
decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
        std::vector<int > &counts);
void
decompose_1d(MPI_Comm &comm, int length, std::vector<int> & offsets,
        std::vector<int> &counts);

int
decompose_1d_local(MPI_Comm &comm, int length);

#endif /* PARALLEL_UTILS_H_ */
