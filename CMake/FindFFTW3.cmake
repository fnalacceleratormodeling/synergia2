SET(FFTW3_INCLUDE_SEARCHPATH
    /usr/local/include
    /usr/local/include/fftw
    /usr/include
    /usr/include/fftw
)
FIND_PATH(FFTW3_INCLUDE_PATH_FOUND fftw3.h ${FFTW3_INCLUDE_SEARCHPATH})
OPTION(FFTW3_INCLUDE_PATH "Include path for FFTW3" ${FFTW3_INCLUDE_PATH_FOUND})

SET(FFTW3_LIB_SEARCHPATH
    /usr/local/lib
    /usr/local/lib/fftw
    /usr/lib
    /usr/lib/fftw
)
FIND_LIBRARY(FFTW3_LIB_FOUND  fftw3 ${FFTW3_LIB_SEARCHPATH})
OPTION(FFTW3_LIB "Library path for FFTW3" ${FFTW3_LIB_FOUND})
FIND_LIBRARY(FFTW3_MPI_LIB_FOUND  fftw3_mpi ${FFTW3_LIB_SEARCHPATH})
OPTION(FFTW3_MPI_LIB "MPI Library path for FFTW3" ${FFTW3_MPI_LIB_FOUND})
