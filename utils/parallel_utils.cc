#include "parallel_utils.h"
#include <cmath>

int
mpi_get_rank(const MPI_Comm &comm)
{
    int error, rank;
    error = MPI_Comm_rank(comm, &rank);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_rank");
    }
    return rank;
}

int
mpi_get_size(const MPI_Comm &comm)
{
    int error, size;
    error = MPI_Comm_size(comm, &size);
    if (error != MPI_SUCCESS) {
        throw std::runtime_error("MPI error in MPI_Comm_size");
    }
    return size;
}

void
decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
        std::vector<int > &counts)
{
    int min_counts = length / processors;
    int remainder = fmod(length, processors);
    int offset = 0;
    for (int i = 0; i < processors; ++i) {
        int count = min_counts;
        if (i >= (processors - remainder)) {
            count += 1;
        }
        offsets.at(i) = offset;
        counts.at(i) = count;
        offset += count;
    }
}

void
decompose_1d(const MPI_Comm &comm, int length, std::vector<int > & offsets,
        std::vector<int > &counts)
{
    int size = mpi_get_size(comm);
    decompose_1d_raw(size, length, offsets, counts);
}

int
decompose_1d_local(const MPI_Comm &comm, int length)
{
    int size = mpi_get_size(comm);
    int rank = mpi_get_rank(comm);
    std::vector<int > offsets(size), counts(size);
    decompose_1d_raw(size, length, offsets, counts);
    return counts[rank];
}
