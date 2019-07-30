#ifndef DISTRIBUTED_FFT2D_FFTW_H_
#define DISTRIBUTED_FFT2D_FFTW_H_

#include <fftw3.h>
#include <fftw3-mpi.h>

#include <vector>
#include <string>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

class Distributed_fft2d
{

private:

    std::array<int, 3> shape;

    MPI_Comm comm;

    fftw_plan plan;
    fftw_plan inv_plan;
    fftw_complex *data;
    fftw_complex *workspace;

    int lower;
    int nx;

public:

    Distributed_fft2d();
    ~Distributed_fft2d();

    int get_lower() const { return lower; }
    int get_upper() const { return lower + nx; }
    std::array<int, 3> const& get_shape() const { return shape; }

    void construct(std::array<int, 3> const& shape, MPI_Comm comm);

    void transform(karray1d_dev & in, karray1d_dev & out);
    void inv_transform(karray1d_dev & in, karray1d_dev & out);
    double get_roundtrip_normalization() const;
};

#endif /* DISTRIBUTED_FFT2D_H_ */
