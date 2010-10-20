#include "parallel_utils.h"
#include <cmath>

void
decompose_1d_raw(int processors, int length, std::vector<int > &offsets,
        std::vector<int > &counts)
{
    int min_counts = length / processors;
    int remainder = static_cast<int > (fmod(length, processors));
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
decompose_1d(Commxx comm, int length, std::vector<int > & offsets, std::vector<
        int > &counts)
{
    int size = comm.get_size();
    decompose_1d_raw(size, length, offsets, counts);
}

int
decompose_1d_local(Commxx comm, int length)
{
    int size = comm.get_size();
    int rank = comm.get_rank();
    std::vector<int > offsets(size), counts(size);
    decompose_1d_raw(size, length, offsets, counts);
    return counts[rank];
}
