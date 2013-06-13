# - FindFFTW3
# Find FFTW3 includes and library
# This module defines:
# FFTW3_INCLUDE_DIR, where to find fftw3.h, etc.
# FFTW3_LIBRARIES, the libraries needed to use fftw3
# FFTW3_LIBRARY_DIRS, the directory containing the fftw3 libraries
# FFTW3_FOUND
#
# also defined are:
# FFTW3_MPI_LIBRARIES, the libraries needed to use parallel fftw3
# FFTW3_MPI_FOUND

# Find the (plain, non-MPI) library
set(FFTW3_NAMES ${FFTW3_NAMES} fftw3)

if(FFTW3_LIBRARY_DIRS)
    find_library(FFTW3_LIBRARIES NAMES ${FFTW3_NAMES} 
        PATHS ${FFTW3_LIBRARY_DIRS}
        NO_DEFAULT_PATH)
endif(FFTW3_LIBRARY_DIRS)

find_library(FFTW3_LIBRARIES NAMES ${FFTW3_NAMES})

if(NOT FFTW3_LIBRARY_DIRS)
    get_filename_component(FFTW3_LIBRARY_DIRS ${FFTW3_LIBRARIES} PATH)
endif(NOT FFTW3_LIBRARY_DIRS)

# Find the include path
if(NOT FFTW3_INCLUDE_DIR)
    get_filename_component(_fftw3_prefix ${FFTW3_LIBRARY_DIRS} PATH)
    find_path(FFTW3_INCLUDE_DIR fftw3.h
        PATHS ${_fftw3_prefix}/include
        NO_DEFAULT_PATH)
endif(NOT FFTW3_INCLUDE_DIR)

find_path(FFTW3_INCLUDE_DIR fftw3.h)

# Find the MPI libraries
set(FFTW3_MPI_NAMES ${FFTW3_MPI_NAMES} fftw3_mpi)

find_library(FFTW3_MPI_LIBRARIES NAMES ${FFTW3_MPI_NAMES}
    PATHS ${FFTW3_LIBRARY_DIRS}
    NO_DEFAULT_PATH)
if(FFTW3_MPI_LIBRARIES)
    set(FFTW3_MPI_FOUND TRUE)
else(FFTW3_MPI_LIBRARIES)
    set(FFTW3_MPI_FOUND FALSE)
endif(FFTW3_MPI_LIBRARIES)

# Find the OpenMP libraries
set(FFTW3_OMP_NAMES ${FFTW3_OMP_NAMES} fftw3_omp)

find_library(FFTW3_OMP_LIBRARIES NAMES ${FFTW3_OMP_NAMES}
    PATHS ${FFTW3_LIBRARY_DIRS}
    NO_DEFAULT_PATH)
if(FFTW3_OMP_LIBRARIES)
    set(FFTW3_OMP_FOUND TRUE)
else(FFTW3_OMP_LIBRARIES)
    set(FFTW3_OMP_FOUND FALSE)
endif(FFTW3_OMP_LIBRARIES)


include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(FFTW3 
    DEFAULT_MSG FFTW3_LIBRARIES FFTW3_INCLUDE_DIR)

