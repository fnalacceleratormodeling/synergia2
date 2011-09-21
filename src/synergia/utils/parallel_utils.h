#ifndef PARALLEL_UTILS_H_
#define PARALLEL_UTILS_H_

#include "mpi.h"
#include <vector>
#include "commxx.h"

/// Perform a one-dimensional decomposition. In cases where the decomposition
/// is uneven, i.e., processors is not an integral multiple of length, the
/// extra items are assigned to the higher-numbered processors.
/// @param processors is the number of processors.
/// @param length is the length of the vector to decompose.
/// @param offsets is a vector of length processors that will be filled with
/// the vector offsets on each processor.
/// @param counts is a vector of length processors that will be filled with
/// the number of counts assigned to each processor.
void
decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
        std::vector<int > &counts);

/// See decompose_1d_raw. The number of processors is extracted from the
/// Commxx object.
void
decompose_1d(Commxx comm, int length, std::vector<int > &offsets,
        std::vector<int > &counts);

/// See decompose_1d_raw. In this simpler version, the number of processors is
/// extracted from the Commxx object and the number of counts assigned to the
/// current processor is returned. The lengths on other processors and offsets are
/// discarded.
int
decompose_1d_local(Commxx comm, int length);

/// Distribute elements across multiple processors. Where decompose assumes that
/// each element will be assigned to a single processor, distribute assumes that
/// elements may be shared among multiple processors.  The return value is a
/// vector of length elements, each member of which is a vector containing the
/// ranks for each element.
/// @param processors is the number of processors.
/// @param elements is the number of elements.
std::vector<std::vector<int > >
distribute_1d_raw(int processors, int elements);

/// See distribute_1d_raw. The number of processors is extracted from the
/// Commxx object.
std::vector<std::vector<int > >
distribute_1d(Commxx comm, int elements);


#endif /* PARALLEL_UTILS_H_ */

