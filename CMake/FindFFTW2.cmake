OPTION(FFTW2_PREFIX "prefix for FFTW2")

IF(FFTW2_PREFIX)
    SET(FFTW2_INCLUDE_SEARCHPATH ${FFTW2_PREFIX}/include)
    SET(FFTW2_LIB_SEARCHPATH ${FFTW2_PREFIX}/lib)
ELSE(FFTW2_PREFIX)
    SET(FFTW2_INCLUDE_SEARCHPATH
        /usr/local/include
        /usr/local/include/fftw
        /usr/include
        /usr/include/fftw
    )
    SET(FFTW2_LIB_SEARCHPATH
        /usr/local/lib
        /usr/local/lib/fftw
        /usr/lib
        /usr/lib/fftw
    )
ENDIF(FFTW2_PREFIX)

FIND_PATH(FFTW2_INCLUDE_PATH fftw.h ${FFTW2_INCLUDE_SEARCHPATH} 
    DOC "Include path for FFTW2")

FIND_LIBRARY(FFTW2_LIB  fftw ${FFTW2_LIB_SEARCHPATH}
    DOC "Library path for FFTW2")
FIND_LIBRARY(FFTW2_MPI_LIB fftw_mpi ${FFTW2_LIB_SEARCHPATH}
    DOC "MPI Library path for FFTW2")
FIND_LIBRARY(FFTW2_R_LIB rfftw ${FFTW2_LIB_SEARCHPATH}
    DOC "Library path for R FFTW2")
FIND_LIBRARY(FFTW2_R_MPI_LIB rfftw_mpi ${FFTW2_LIB_SEARCHPATH}
    DOC "MPI Library path for R FFTW2")

#~ MESSAGE("jfa debug: FFTW2_INCLUDE_PATH=${FFTW2_INCLUDE_PATH}")
#~ MESSAGE("jfa debug: FFTW2_LIB=${FFTW2_LIB}")