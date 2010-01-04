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
message("FFTW3_ROOT = '${FFTW3_ROOT}'")
message("FFTW3_INCLUDE_SEARCHPATH = '${FFTW3_INCLUDE_SEARCHPATH}'")
FIND_PATH(FFTW3_INCLUDE_PATH_FOUND fftw3.h ${FFTW3_INCLUDE_SEARCHPATH})
#~ OPTION(FFTW3_INCLUDE_PATH "Include path for FFTW3" ${FFTW3_INCLUDE_PATH_FOUND})
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
FIND_LIBRARY(FFTW3_LIB_FOUND  fftw3 ${FFTW3_LIB_SEARCHPATH})
#~ OPTION(FFTW3_LIB "Library path for FFTW3" ${FFTW3_LIB_FOUND})
set(FFTW3_LIB ${FFTW3_LIB_FOUND})
#~ FIND_LIBRARY(FFTW3_MPI_LIB_FOUND  fftw3_mpi ${FFTW3_LIB_SEARCHPATH})
#~ OPTION(FFTW3_MPI_LIB "MPI Library path for FFTW3" ${FFTW3_MPI_LIB_FOUND})

#~ MESSAGE("jfa debug: FFTW3_INCLUDE_PATH=${FFTW3_INCLUDE_PATH}")
#~ MESSAGE("jfa debug: FFTW3_INCLUDE_PATH_FOUND=${FFTW3_INCLUDE_PATH_FOUND}")
#~ MESSAGE("jfa debug: FFTW3_LIB=${FFTW3_LIB}")
#~ MESSAGE("jfa debug: FFTW3_MPI_LIB=${FFTW3_MPI_LIB}")