project(SYNERGIA2)
cmake_minimum_required(VERSION 2.8)
enable_testing()
set(CMAKE_MODULE_PATH "${SYNERGIA2_SOURCE_DIR}/CMake")

# MPI
FIND_PACKAGE(MPI REQUIRED)
INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_CXX_COMPILER(${MPI_COMPILER} "MPI C++ Compiler")

# python
find_package(PythonInterp REQUIRED)
find_package(PythonLibs REQUIRED)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})

# boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD ON)
FIND_PACKAGE(Boost COMPONENTS python regex unit_test_framework serialization REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# hdf5
find_package(HDF5 REQUIRED)

# chef
find_package(CHEF REQUIRED)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})

# numpy
find_package(NUMPY REQUIRED)
include_directories(${NUMPY_INCLUDE_DIR})

# mpi4py
find_package(MPI4PY REQUIRED)
include_directories(${MPI4PY_INCLUDE_DIR})

# misc
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindGSL.cmake)

# fftw3 (and, possibly, fftw2)
find_package(FFTW3 REQUIRED)
if(FFTW3_MPI_FOUND)
    include_directories(${FFTW3_INCLUDE_DIR})
    set(PARALLEL_FFTW_LIBRARIES ${FFTW3_MPI_LIBRARIES} ${FFTW3_LIBRARIES})
    if(NOT PARALLEL_FFTW_FOUND)
        message("-- Using FFTW3 for parallel FFTs")
        set(PARALLEL_FFTW_FOUND TRUE CACHE BOOL "found parallel fftw" FORCE)
    endif(NOT PARALLEL_FFTW_FOUND)
else(FFTW3_MPI_FOUND)
    message("-- FFTW3 MPI libraries not found, trying FFTW2")
    find_package(FFTW2)
    if(NOT FFTW2_MPI_FOUND)
        message(FATAL_ERROR "Synergia requires either FFTW2 or FFT3 compiled with MPI support.")
    endif(NOT FFTW2_MPI_FOUND)
    include_directories(${FFTW2_INCLUDE_DIR})
    set(PARALLEL_FFTW_LIBRARIES ${FFTW2_MPI_LIBRARIES} ${FFTW2_LIBRARIES})
    add_definitions(-DUSE_FFTW2)
    if(NOT PARALLEL_FFTW_FOUND)
        message("-- Using FFTW2 for parallel FFTs")
        set(PARALLEL_FFTW_FOUND TRUE CACHE BOOL "found parallel fftw" FORCE)
    endif(NOT PARALLEL_FFTW_FOUND)
endif(FFTW3_MPI_FOUND)

# for eclipse only
IF(CMAKE_COMPILER_IS_GNUCC)
  SET(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCC)
IF(CMAKE_COMPILER_IS_GNUCXX)
  SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fmessage-length=0")
ENDIF(CMAKE_COMPILER_IS_GNUCXX)

# External libraries generate many warnings
#add_definitions(-Wall -pedantic -Wno-long-long)

set(CHEF_LIBRARY_DIRS_EXPANDED "")
foreach(dir ${CHEF_LIBRARY_DIRS})
    set(CHEF_LIBRARY_DIRS_EXPANDED "${CHEF_LIBRARY_DIRS_EXPANDED} \"${dir}\"")
endforeach(dir)

configure_file("${SYNERGIA2_SOURCE_DIR}/synergia-local.in"
    "${SYNERGIA2_BINARY_DIR}/synergia-local" IMMEDIATE)

configure_file("${SYNERGIA2_SOURCE_DIR}/src/local_paths.py.in"
    "${SYNERGIA2_BINARY_DIR}/src/local_paths.py" IMMEDIATE)

include_directories(${SYNERGIA2_SOURCE_DIR}/src)

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(docs)
