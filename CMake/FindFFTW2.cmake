# - FindFFTW2
# Find FFTW2 includes and library
# This module defines:
# FFTW2_INCLUDE_DIR, where to find fftw2.h, etc.
# FFTW2_LIBRARIES, the libraries needed to use fftw2
# FFTW2_LIBRARY_DIRS, the directory containing the fftw2 libraries
# FFTW2_FOUND
#
# also defined are:
# FFTW2_MPI_LIBRARIES, the libraries needed to use parallel fftw2
# FFTW2_MPI_FOUND

# Find the (plain, non-MPI) libraries
set(FFTW2_NAME ${FFTW2_NAME} fftw)
set(FFTW2_RNAME ${FFTW2_RNAME} rfftw)

if(FFTW2_LIBRARY_DIRS)
    find_library(FFTW2_LIBRARY NAMES ${FFTW2_NAME} 
        PATHS ${FFTW2_LIBRARY_DIRS}
        NO_DEFAULT_PATH)
    find_library(FFTW2_RLIBRARY NAMES ${FFTW2_RNAME} 
        PATHS ${FFTW2_LIBRARY_DIRS}
        NO_DEFAULT_PATH)
endif(FFTW2_LIBRARY_DIRS)

find_library(FFTW2_LIBRARY NAMES ${FFTW2_NAME})
find_library(FFTW2_RLIBRARY NAMES ${FFTW2_RNAME})
set(FFTW2_LIBRARIES ${FFTW2_LIBRARY} ${FFTW2_RLIBRARY})

if(NOT FFTW2_LIBRARY_DIRS)
    get_filename_component(FFTW2_LIBRARY_DIRS ${FFTW2_LIBRARY} PATH)
endif(NOT FFTW2_LIBRARY_DIRS)

message("FFTW2_LIBRARY_DIRS=${FFTW2_LIBRARY_DIRS}")
# Find the include path
if(NOT FFTW2_INCLUDE_DIR AND FFTW2_LIBRARY_DIRS)
    get_filename_component(_fftw2_prefix ${FFTW2_LIBRARY_DIRS} PATH)
    find_path(FFTW2_INCLUDE_DIR fftw.h
        PATHS ${_fftw2_prefix}/include
        NO_DEFAULT_PATH)
endif(NOT FFTW2_INCLUDE_DIR AND FFTW2_LIBRARY_DIRS)

find_path(FFTW2_INCLUDE_DIR fftw.h)

# Find the MPI libraries
set(FFTW2_MPI_NAME ${FFTW2_MPI_NAME} fftw_mpi)
set(FFTW2_MPI_RNAME ${FFTW2_MPI_RNAME} rfftw_mpi)

find_library(FFTW2_MPI_LIBRARY NAMES ${FFTW2_MPI_NAME}
    PATHS ${FFTW2_LIBRARY_DIRS}
    NO_DEFAULT_PATH)
find_library(FFTW2_MPI_RLIBRARY NAMES ${FFTW2_MPI_RNAME}
    PATHS ${FFTW2_LIBRARY_DIRS}
    NO_DEFAULT_PATH)
set(FFTW2_MPI_LIBRARIES ${FFTW2_MPI_LIBRARY} ${FFTW2_MPI_RLIBRARY})
if(FFTW2_MPI_LIBRARIES)
    set(FFTW2_MPI_FOUND TRUE)
else(FFTW2_MPI_LIBRARIES)
    set(FFTW2_MPI_FOUND FALSE)
endif(FFTW2_MPI_LIBRARIES)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW2 
    DEFAULT_MSG FFTW2_LIBRARIES FFTW2_INCLUDE_DIR)

