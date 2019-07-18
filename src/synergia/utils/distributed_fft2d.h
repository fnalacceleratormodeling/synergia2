#ifndef DISTRIBUTED_FFT2D
#define DISTRIBUTED_FFT2D

#include <Kokkos_Core.hpp>

#ifdef KOKKOS_ENABLE_CUDA

  #include "synergia/utils/distributed_fft2d_cuda.h"

#else

  #include "synergia/utils/distributed_fft2d_fftw.h"

#endif


#endif
