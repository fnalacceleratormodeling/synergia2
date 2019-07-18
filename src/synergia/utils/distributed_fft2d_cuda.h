#ifndef DISTRIBUTED_FFT2D_H_
#define DISTRIBUTED_FFT2D_H_

#include <vector>
#include <string>

#include <cufft.h>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

class Distributed_fft2d
{

private:

    std::array<int, 3> shape;
    MPI_Comm comm;

    cufftHandle plan;

    int lower, upper;
    std::vector<int> uppers, lengths, lengths_1d;

public:

    Distributed_fft2d(std::array<int, 3> const& shape, MPI_Comm comm);
    ~Distributed_fft2d();

    int get_lower() const;
    int get_upper() const;

    std::vector<int>   const& get_uppers();
    std::vector<int>   const& get_lengths();
    std::vector<int>   const& get_lengths_1d();
    std::array<int, 3> const& get_shape() const;

    void transform(karray1d_dev & in, karray1d_dev & out);
    void inv_transform(karray1d_dev & in, karray1d_dev & out);
    double get_roundtrip_normalization() const;
};

#endif /* DISTRIBUTED_FFT2D_H_ */
