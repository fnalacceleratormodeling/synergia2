#ifndef DISTRIBUTED_FFT3D_CUDA_H_
#define DISTRIBUTED_FFT3D_CUDA_H_

#include <vector>
#include <string>

#include <cufft.h>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

class Distributed_fft3d
{

private:

    std::array<int, 3> shape;
    Commxx comm;

    cufftHandle plan;
    cufftHandle invplan;

    int lower;
    int nz;

public:

    static int get_padded_shape_real(int s)
    { return 2*(s/2+1); }

    static int get_padded_shape_cplx(int s)
    { return s/2+1; }

    Distributed_fft3d();
    ~Distributed_fft3d();


    int get_lower() const { return lower; }
    int get_upper() const { return lower + nz; }

    int padded_nx_real() const { return get_padded_shape_real(shape[0]); }
    int padded_nx_cplx() const { return get_padded_shape_cplx(shape[0]); }

    std::array<int, 3> const& get_shape() const { return shape; }
    Commxx const& get_comm() const { return comm; }

    void construct(std::array<int, 3> const& shape, Commxx const& comm);

    void transform(karray1d_dev& in, karray1d_dev& out);
    void inv_transform(karray1d_dev& in, karray1d_dev& out);
    double get_roundtrip_normalization() const;
};

#endif /* DISTRIBUTED_FFT2D_H_ */
