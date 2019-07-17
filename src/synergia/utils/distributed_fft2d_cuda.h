#ifndef DISTRIBUTED_FFT2D_H_
#define DISTRIBUTED_FFT2D_H_

#include <fftw3.h>
#include <fftw3-mpi.h>

#include <vector>
#include <string>

#include "synergia/utils/multi_array_typedefs.h"
#include "synergia/utils/commxx.h"

template<typename FFT>
class Distributed_fft2d_base
{

private:

    std::array<int, 3> shape;

    FFT fft;

#if 0
    fftw_plan plan, inv_plan;
    fftw_complex *data;
    fftw_complex *workspace;
#endif

    int lower, upper;
    std::vector<int> uppers, lengths, lengths_1d;

    int  local_size_real;
    bool have_local_data;

    void calculate_uppers_lengths();

public:

    Distributed_fft2d(std::array<int, 3> const& shape);
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

#ifdef KOKKOS_HAVE_CUDA

  // cuda implementation of FFT
  #include "fft2d_impl_cuda.h"
  typedef Distributed_fft2d_base<fft2d_impl_cuda> Distributed_fft2d;

#else

  // FFTW
  #include "fft2d_impl_fftw.h"
  typedef Distributed_fft2d_base<fft2d_impl_fftw> Distributed_fft2d;

#endif



#endif /* DISTRIBUTED_FFT2D_H_ */
