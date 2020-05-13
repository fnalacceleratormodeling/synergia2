#ifndef DISTRIBUTED_FFT3D
#define DISTRIBUTED_FFT3D

#include <Kokkos_Core.hpp>

#ifdef Kokkos_ENABLE_CUDA

  #include "synergia/utils/distributed_fft3d_cuda.h"

#else

  #include "synergia/utils/distributed_fft3d_fftw.h"

#endif


#endif
