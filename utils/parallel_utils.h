#ifndef PARALLEL_UTILS_H_
#define PARALLEL_UTILS_H_

#include "mpi.h"
#include <vector>
#include "mpi_comm_wrap.h"

void
decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
        std::vector<int > &counts);

void
decompose_1d(MPI_comm_wrap comm, int length, std::vector<int > &offsets,
        std::vector<int > &counts);

int
decompose_1d_local(MPI_comm_wrap comm, int length);

#endif /* PARALLEL_UTILS_H_ */
