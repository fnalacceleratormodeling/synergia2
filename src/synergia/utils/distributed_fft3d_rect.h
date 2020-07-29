#ifndef DISTRIBUTED_FFT3D_H
#define DISTRIBUTED_FFT3D_H

// discrete FFT on 3D grid. 
// 2D DST-II on x-y plane and DFT along the z-axis

#include <Kokkos_Core.hpp>

#ifdef Kokkos_ENABLE_CUDA

  #include "synergia/utils/distributed_fft3d_rect_cuda.h"

#else

  #include "synergia/utils/distributed_fft3d_rect_fftw.h"

#endif


#endif
