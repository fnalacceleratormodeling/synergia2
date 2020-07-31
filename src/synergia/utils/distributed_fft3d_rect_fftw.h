#ifndef DISTRIBUTED_FFT3D_RECT_FFTW_H_
#define DISTRIBUTED_FFT3D_RECT_FFTW_H_

#include <fftw3.h>
#include <fftw3-mpi.h>

#include <array>
#include <string>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

class Distributed_fft3d_rect
{

private:

    std::array<int, 3> shape;

    MPI_Comm comm;

    fftw_plan plan_xy, plan_z;
    fftw_plan inv_plan_xy, inv_plan_z;

    double *data;
    fftw_complex *workspace;

    int lower;
    int nx;

public:

    static int get_padded_shape_real(int s)
    { return 2*(s/2+1); }

    static int get_padded_shape_cplx(int s)
    { return s/2+1; }

    Distributed_fft3d_rect();
    ~Distributed_fft3d_rect();

    int get_lower() const { return lower; }
    int get_upper() const { return lower + nx; }

    int padded_nz_real() const { return get_padded_shape_real(shape[2]); }
    int padded_nz_cplx() const { return get_padded_shape_cplx(shape[2]); }

    std::array<int, 3> const& get_shape() const { return shape; }

    void construct(std::array<int, 3> const& shape, MPI_Comm comm);

    void transform(karray1d_dev & in, karray1d_dev & out);
    void inv_transform(karray1d_dev & in, karray1d_dev & out);
    double get_roundtrip_normalization() const;
};

#endif /* DISTRIBUTED_FFT3D_H_ */
