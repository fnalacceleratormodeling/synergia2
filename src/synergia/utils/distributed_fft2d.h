#ifndef DISTRIBUTED_FFT2D
#define DISTRIBUTED_FFT2D

#if KOKKOS_HAVE_CUDA

  #include "synergia/utils/distributed_fft2d_cuda.h"

#else

  #include "synergia/utils/distributed_fft2d_fftw.h"

#endif


#endif
