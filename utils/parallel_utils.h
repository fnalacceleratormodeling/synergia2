#ifndef PARALLEL_UTILS_H_
#define PARALLEL_UTILS_H_

#include "mpi.h"
#include <vector>
#include "commxx.h"

void
decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
        std::vector<int > &counts);

void
decompose_1d(Commxx comm, int length, std::vector<int > &offsets,
        std::vector<int > &counts);

int
decompose_1d_local(Commxx comm, int length);

#endif /* PARALLEL_UTILS_H_ */
