project(SYNERGIA2)
set(SYNERGIA_MAJOR_VERSION 2)
set(SYNERGIA_MINOR_VERSION 0)
set(SYNERGIA_SUBMINOR_VERSION 90)

cmake_minimum_required(VERSION 2.8)
enable_testing()
set(CMAKE_MODULE_PATH "${SYNERGIA2_SOURCE_DIR}/CMake")

# Set defaults for building python bindings (build bindings by default).
if (NOT DEFINED BUILD_PYTHON_BINDINGS)
  if (NOT $ENV{BUILD_PYTHON_BINDINGS} STREQUAL "")
    set (BUILD_PYTHON_BINDINGS $ENV{BUILD_PYTHON_BINDINGS})
  else ()
    set (BUILD_PYTHON_BINDINGS TRUE)
  endif ()
endif ()

if (NOT BUILD_PYTHON_BINDINGS)
  message (STATUS "Not building python bindings")
endif ()

# MPI
FIND_PACKAGE(MPI REQUIRED)
INCLUDE(CMakeForceCompiler)
CMAKE_FORCE_CXX_COMPILER(${MPI_COMPILER} "MPI C++ Compiler")

# python
if (BUILD_PYTHON_BINDINGS)
    find_package(PythonInterp REQUIRED)
    find_package(PythonLibs REQUIRED)
    INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
endif()

# boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD ON)
FIND_PACKAGE(Boost COMPONENTS python regex unit_test_framework serialization REQUIRED)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# hdf5
find_package(HDF5 REQUIRED CXX)
INCLUDE_DIRECTORIES(${HDF5_INCLUDE_DIRS})

# chef
find_package(CHEF REQUIRED)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})

# numpy
find_package(NUMPY REQUIRED)
include_directories(${NUMPY_INCLUDE_DIR})

if (BUILD_PYTHON_BINDINGS)
    # mpi4py
    find_package(MPI4PY REQUIRED)
    include_directories(${MPI4PY_INCLUDE_DIR})
endif()

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

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DBOOST_DISABLE_ASSERTS")
if(USE_SIMPLE_TIMER)
    message("-- Simple timer profiling enabled")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_SIMPLE_TIMER")
endif(USE_SIMPLE_TIMER)

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

include_directories(BEFORE ${SYNERGIA2_SOURCE_DIR}/src)

# Are we building shared or static (shared is the default)?
if (NOT DEFINED BUILD_SHARED_LIBS)
  if (NOT $ENV{BUILD_SHARED_LIBS} STREQUAL "")
    set (BUILD_SHARED_LIBS $ENV{BUILD_SHARED_LIBS})
  else ()
    set (BUILD_SHARED_LIBS TRUE)
  endif ()
endif ()

if (BUILD_SHARED_LIBS)
  message (STATUS "Building shared libraries")
else ()
    # This prevents multiple-definition errors.
    set(CMAKE_EXE_LINKER_FLAGS
         "${CMAKE_EXE_LINKER_FLAGS} -Wl,--allow-multiple-definition")
  message (STATUS "Building static libraries")
endif ()

# We should probably warn someone trying to build static libraries along with
#  the python bindings since it won't work with python.
# The python bindings will still be build shared.
if ( (NOT BUILD_SHARED_LIBS) AND BUILD_PYTHON_BINDINGS )
  message(WARNING "Static library build selected, ignored for python bindings")
endif ()

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(docs)
add_subdirectory(synergia-script-templates)

# dummy comment -- please delete me
