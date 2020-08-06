
project(SYNERGIA2)

cmake_minimum_required(VERSION 3.10)
enable_testing()

# project specific cmake modules
set(CMAKE_MODULE_PATH "${SYNERGIA2_SOURCE_DIR}/CMake")

# C++ 14 standard
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# versions
set(SYNERGIA2_VERSION_MAJOR 2018) # year
set(SYNERGIA2_VERSION_MINOR 02)   # month
set(SYNERGIA2_VERSION_PATCH 20)   # day
set(SYNERGIA2_VERSION_TWEAK 00)   # sub-day counter
set(SYNERGIA2_VERSION "${SYNERGIA2_VERSION_MAJOR}.${SYNERGIA2_VERSION_MINOR}.${SYNERGIA2_VERSION_PATCH}.${SYNERGIA2_VERSION_TWEAK}")

set(SYNERGIA_MAJOR_VERSION 3)
set(SYNERGIA_MINOR_VERSION 1)
set(SYNERGIA_SUBMINOR_VERSION 00)


# needed to avoid warnings about rpath policy on Mac OSX
# no known difference between "OLD" and "NEW" policies for us
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

# misc
include(CopyFile)
include(MPITest)
include(AddPythonExtension)
include(SmartIncludeDirectories)

# include dirs
set(extra_includes)

# options
option(SIMPLE_TIMER "use the simple timer 2 profiling" Off)
option(BUILD_PYTHON_BINDINGS "build the Python bindings of Synergia" On)
option(ALLOW_PADDING "add paddings to the particle array for memory alignment" On)

# only build the shared lib
option(BUILD_SHARED_LIBS "Build shared libraries" TRUE)

# find MPI -- do not build the C++ bindings
set(MPI_CXX_SKIP_MPICXX TRUE CACHE BOOL "disable the MPI2-C++ bindings")
find_package(MPI REQUIRED COMPONENTS C CXX) # CXX means we need to call C-APIs in C++ code
list(APPEND extra_includes ${MPI_C_INCLUDE_DIRS})
set_property(DIRECTORY PROPERTY COMPILE_DEFINITIONS ${MPI_CXX_COMPILE_DEFINITIONS})
execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE MPI_VERSION)

# mpi oversubscribe
if(${MPI_VERSION} MATCHES "OpenRTE.* 3")
    set(MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
else ()
    set(MPI_OVERSUBSCRIBE_FLAG "")
endif()
set(MPIEXEC_NUMPROC_FLAG ${MPI_OVERSUBSCRIBE_FLAG} ${MPIEXEC_NUMPROC_FLAG})

# OpenMP
if(${Kokkos_ENABLE_OPENMP})
    FIND_PACKAGE(OpenMP)
    if (OPENMP_FOUND)
        set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
        set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
    endif ()
endif()

# hdf5
set(HDF5_PREFER_PARALLEL true)
find_package(HDF5)
list(APPEND extra_includes ${HDF5_INCLUDE_DIRS})

if(HDF5_IS_PARALLEL)
    message(STATUS "Using parallel hdf5")
    add_definitions(-DUSE_PARALLEL_HDF5)
else(HDF5_IS_PARALLEL)
    message(STATUS "Using serial hdf5")
endif()

# GSL
find_package(GSL REQUIRED)
list(APPEND extra_includes ${GSL_INCLUDE_DIR})

# Eigen3
find_package(Eigen3 REQUIRED)
list(APPEND extra_includes ${EIGEN3_INCLUDE_DIR})

# fftw3
if(${Kokkos_ENABLE_OPENMP})
    find_package(FFTW3 REQUIRED)
    list(APPEND extra_includes ${FFTW3_INCLUDE_DIR})

    if(NOT FFTW3_MPI_FOUND)
        message(FATAL_ERROR "Missing FFTW3-mpi")
    endif()

    if (OPENMP_FOUND AND NOT FFTW3_OMP_FOUND)
        message(FATAL_ERROR "FFTW3 openmp not found")
    endif()

    if(OPENMP_FOUND)
        set(PARALLEL_FFTW_LIBRARIES ${FFTW3_MPI_LIBRARIES} ${FFTW3_OMP_LIBRARIES} ${FFTW3_LIBRARIES})
    else(OPENMP_FOUND)
        set(PARALLEL_FFTW_LIBRARIES ${FFTW3_MPI_LIBRARIES} ${FFTW3_LIBRARIES})
    endif(OPENMP_FOUND)
endif()

# warnings
string(APPEND ${CMAKE_CXX_FLAGS} "-Wall -Werror")

# enable fast-math, but disable the finite-math-only to avoid the std::isfinite() woes
string(APPEND ${CMAKE_CXX_FLAGS} "-ffast-math -fno-finite-math-only")

# for eclipse only (no line wrapping for error messages)
IF(CMAKE_COMPILER_IS_GNUCXX)
    string(APPEND ${CMAKE_CXX_FLAGS} "-fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# python
if (BUILD_PYTHON_BINDINGS)
    message (STATUS "Building python bindings")
    add_subdirectory(src/synergia/utils/pybind11)

    #find_package(ConsistentPython)
    #list(APPEND extra_includes ${PYTHON_INCLUDE_DIR})
    #set(IPYTHON_EXECUTABLE "ipython" CACHE PATH "path to optional ipython executable")

    # mpi4py
    #find_package(MPI4PY REQUIRED)
    #list(APPEND extra_includes ${MPI4PY_INCLUDE_DIR})

else(BUILD_PYTHON_BINDINGS)
    message (STATUS "Do not build python bindings")
endif()


# simple timer
if(SIMPLE_TIMER)
    message(STATUS "Simple timer 2 profiling enabled")
    add_definitions(-DSIMPLE_TIMER)
    if (SIMPLE_TIMER_BARRIER)
        message(STATUS "Simple timer 2 MPI_Barrier enabled")
        add_definitions(-DSIMPLE_TIMER_BARRIER)
    endif()
endif()

# allow padding
if(ALLOW_PADDING)
    message(STATUS "Particle array padding enabled")
else()
    message(STATUS "Particle array padding disabled")
    add_definitions(-DNO_PADDING)
endif()

# Additional CXXFLAGS
if (EXTRA_CXX_FLAGS)
    string(APPEND CMAKE_CXX_FLAGS ${EXTRA_CXX_FLAGS})
endif()

# External libraries generate many warnings
#add_definitions(-Wall -pedantic -Wno-long-long)

configure_file("${SYNERGIA2_SOURCE_DIR}/synergia-local.in"
    "${SYNERGIA2_BINARY_DIR}/synergia-local" IMMEDIATE)

configure_file("${SYNERGIA2_SOURCE_DIR}/src/local_paths.py.in"
    "${SYNERGIA2_BINARY_DIR}/src/local_paths.py" IMMEDIATE)

include_directories(BEFORE ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils)
include_directories(BEFORE ${SYNERGIA2_SOURCE_DIR}/src)


if(NOT TARGET test-exes)
  add_custom_target(test-exes)
endif()


function(add_test_executable name)
#  add_executable(${name} EXCLUDE_FROM_ALL ${ARGN})
  add_executable(${name} ${ARGN})
  add_dependencies(test-exes ${name})
endfunction()

add_custom_target(check
  COMMAND ${CMAKE_CTEST_COMMAND}
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
  COMMENT "Running ctest")
add_dependencies(check test-exes)

set(INCLUDE_INSTALL_DIR include/ CACHE PATH "include install directory")
set(LIB_INSTALL_DIR lib/ CACHE PATH "library install directory")
set(BIN_INSTALL_DIR bin/ CACHE PATH "executable install directory")
set(PYTHON_INSTALL_DIR "lib/python${PYTHON_VERSION_MAJOR}.${PYTHON_VERSION_MINOR}/site-packages/"
    CACHE PATH "python install directory")

set(SYNERGIA_EXTRA_LIBRARY_PATHS "" CACHE STRING "additional library paths for synergia executable")
set(SYNERGIA_EXTRA_PYTHON_PATHS "" CACHE STRING "additional python module paths for synergia executable")

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

# include paths
smart_include_directories("${extra_includes}")

# kokkos
add_subdirectory(src/synergia/utils/kokkos)
include_directories(src/synergia/utils/kokkos/core/src)
include_directories(${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos)

if(${Kokkos_ENABLE_CUDA})
    add_definitions(-DKokkos_ENABLE_CUDA)
endif()

if(${Kokkos_ENABLE_OPENMP})
    add_definitions(-DKokkos_ENABLE_OPENMP)
endif()

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(archived-applications)
#add_subdirectory(docs)
add_subdirectory(synergia-script-templates)

