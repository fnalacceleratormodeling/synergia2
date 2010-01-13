project(SYNERGIA2)
cmake_minimum_required(VERSION 2.6)
enable_testing()

# MPI
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindMPIsynergia.cmake)
FIND_PACKAGE(MPI REQUIRED)
INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_CXX_COMPILER(${MPI_COMPILER} "MPI C++ Compiler")

# python
INCLUDE(FindPythonInterp)
INCLUDE(FindPythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
LINK_DIRECTORIES(${PYTHON_LIBRARIES})

# boost
SET(Boost_USE_STATIC_LIBS OFF CACHE BOOL "Do not modify")
SET(Boost_USE_MULTITHREADED ON CACHE BOOL "Do not modify")
SET(Boost_ADDITIONAL_VERSIONS 1.38 1.39 1.40 1.41 1.42 1.43)
FIND_PACKAGE(Boost COMPONENTS python regex unit_test_framework)
message("Boost_INCLUDE_DIRS = '${Boost_INCLUDE_DIRS}'")
message("Boost_LIBRARIES = '${Boost_LIBRARIES}'")
message("Boost_LIBRARY_DIRS = '${Boost_LIBRARY_DIRS}'")
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# chef
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindCHEF.cmake)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
if(NOT CHEF_FOUND)
    message(FATAL_ERROR "CHEF installation not found. Please make sure chef-config.sh is in your path.")
endif(NOT CHEF_FOUND)
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})

# numpy
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindNumpy.cmake)
include_directories(${NUMPY_INCLUDE_DIR})

# misc
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindGSL.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindFFTW2.cmake)

# fftw3
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindFFTW3.cmake)
INCLUDE_DIRECTORIES(${FFTW3_INCLUDE_PATH})

configure_file("${SYNERGIA2_SOURCE_DIR}/synergia-local.in" 
    "${SYNERGIA2_BINARY_DIR}/synergia-local" IMMEDIATE)

include_directories(${SYNERGIA2_SOURCE_DIR})

add_subdirectory(array_nd)
add_subdirectory(synergia)
add_subdirectory(s2_fish)
add_subdirectory(impact)
add_subdirectory(frontend)

add_subdirectory(octapy)

# jfa: hardcoded!!!!
include_directories("/usr/include/glib-2.0")
include_directories("/usr/lib/glib-2.0/include")

add_subdirectory(examples/mu2e)
add_subdirectory(tests/dependencies)
