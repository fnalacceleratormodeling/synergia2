SET(FFTW3_INCLUDE_SEARCHPATH
    /usr/local/include
    /usr/local/include/fftw
    /usr/include
    /usr/include/fftw
)
option(FFTW3_ROOT "Install prefix for fftw3")
if (FFTW3_ROOT)
    set(FFTW3_INCLUDE_SEARCHPATH ${FFTW3_ROOT}/include)
endif (FFTW3_ROOT)
FIND_PATH(FFTW3_INCLUDE_PATH_FOUND fftw3.h HINTS ${FFTW3_INCLUDE_SEARCHPATH})
set(FFTW3_INCLUDE_PATH ${FFTW3_INCLUDE_PATH_FOUND})

SET(FFTW3_LIB_SEARCHPATH
    /usr/local/lib
    /usr/local/lib/fftw
    /usr/lib
    /usr/lib/fftw
)
if (FFTW3_ROOT)
    set(FFTW3_LIB_SEARCHPATH ${FFTW3_ROOT}/lib)
endif (FFTW3_ROOT)
FIND_LIBRARY(FFTW3_LIB_FOUND fftw3 HINTS ${FFTW3_LIB_SEARCHPATH})
set(FFTW3_LIB ${FFTW3_LIB_FOUND})
FIND_LIBRARY(FFTW3_MPI_LIB_FOUND fftw3_mpi HINTS ${FFTW3_LIB_SEARCHPATH})
set(FFTW3_MPI_LIB  ${FFTW3_MPI_LIB_FOUND})
