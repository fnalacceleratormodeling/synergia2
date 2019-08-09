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

    int lower;
    int nx;

public:

    Distributed_fft2d();
    ~Distributed_fft2d();

    void construct(std::array<int, 3> const& shape, MPI_Comm comm);

    int get_lower() const { return lower; }
    int get_upper() const { return lower + nx; }
    std::array<int, 3> const& get_shape() const { return shape; }

    void transform(karray1d_dev & in, karray1d_dev & out);
    void inv_transform(karray1d_dev & in, karray1d_dev & out);
    double get_roundtrip_normalization() const;
};

#endif /* DISTRIBUTED_FFT2D_H_ */
