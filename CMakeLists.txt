project(SYNERGIA2)
if (CMAKE_VERSION VERSION_GREATER 3.0.0)
    set(CMAKE_CXX_STANDARD 14)
else ()
    add_compile_options(-std=c++14)
endif ()

set(SYNERGIA2_VERSION_MAJOR 2018) # year
set(SYNERGIA2_VERSION_MINOR 02)   # month
set(SYNERGIA2_VERSION_PATCH 20)   # day
set(SYNERGIA2_VERSION_TWEAK 00)   # sub-day counter
set(SYNERGIA2_VERSION
    "${SYNERGIA2_VERSION_MAJOR}.${SYNERGIA2_VERSION_MINOR}.${SYNERGIA2_VERSION_PATCH}.${SYNERGIA2_VERSION_TWEAK}")
set(SYNERGIA_MAJOR_VERSION 3)
set(SYNERGIA_MINOR_VERSION 1)
set(SYNERGIA_SUBMINOR_VERSION 00)


cmake_minimum_required(VERSION 2.8)
enable_testing()
set(CMAKE_MODULE_PATH "${SYNERGIA2_SOURCE_DIR}/CMake")

# needed to avoid warnings about rpath policy on Mac OSX
# no known difference between "OLD" and "NEW" policies for us
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

include(SmartIncludeDirectories)
set(extra_includes)

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
FIND_PACKAGE(MPI REQUIRED CXX)
execute_process(COMMAND ${MPIEXEC} --version OUTPUT_VARIABLE MPI_VERSION)
if(${MPI_VERSION} MATCHES "OpenRTE.* 3")
    set(MPI_OVERSUBSCRIBE_FLAG "--oversubscribe")
else ()
    set(MPI_OVERSUBSCRIBE_FLAG "")
endif()
set(MPIEXEC_NUMPROC_FLAG ${MPI_OVERSUBSCRIBE_FLAG} ${MPIEXEC_NUMPROC_FLAG})
if (CMAKE_VERSION VERSION_LESS 3.0.0)
    INCLUDE(CMakeForceCompiler)
    CMAKE_FORCE_CXX_COMPILER(${MPI_COMPILER} "MPI C++ Compiler")
else()
    # jfa: this is really ugly, but provides compatibility 
    # with both CMake 2 and 3
    include_directories(${MPI_CXX_INCLUDE_PATH})
    link_libraries(${MPI_CXX_LIBRARIES})
endif ()

# OpenMP
FIND_PACKAGE(OpenMP)
if (OPENMP_FOUND)
    set (CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
    set (CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
endif ()

# python
if (BUILD_PYTHON_BINDINGS)
    find_package(ConsistentPython)
    list(APPEND extra_includes ${PYTHON_INCLUDE_DIR})
    set(IPYTHON_EXECUTABLE "ipython" CACHE PATH "path to optional ipython executable")
endif()

# boost
SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD ON)
SET(Boost_NO_BOOST_CMAKE ON)
FIND_PACKAGE(Boost COMPONENTS python regex unit_test_framework serialization
    filesystem system iostreams REQUIRED)
list(APPEND extra_includes ${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

# hdf5
find_package(HDF5)
list(APPEND extra_includes ${HDF5_INCLUDE_DIRS})

# chef
find_package(CHEF REQUIRED)
list(APPEND extra_includes ${CHEF_INCLUDE_DIR})
LINK_DIRECTORIES(${CHEF_LIB_DIR})

# numpy
find_package(NUMPY REQUIRED)
include_directories(${NUMPY_INCLUDE_DIR})

if (BUILD_PYTHON_BINDINGS)
    # mpi4py
    find_package(MPI4PY REQUIRED)
    list(APPEND extra_includes ${MPI4PY_INCLUDE_DIR})
endif()

# misc
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindGSL.cmake)
list(APPEND extra_includes ${GSL_INCLUDE_DIR})

# Eigen3
find_package(Eigen3 REQUIRED)
list(APPEND extra_includes ${EIGEN3_INCLUDE_DIR})

# fftw3 (and, possibly, fftw2)
find_package(FFTW3 REQUIRED)
if(FFTW3_MPI_FOUND)
    list(APPEND extra_includes ${FFTW3_INCLUDE_DIR})
    if(OPENMP_FOUND)
        set(PARALLEL_FFTW_LIBRARIES ${FFTW3_MPI_LIBRARIES} ${FFTW3_OMP_LIBRARIES} ${FFTW3_LIBRARIES})
    else(OPENMP_FOUND)
    set(PARALLEL_FFTW_LIBRARIES ${FFTW3_MPI_LIBRARIES} ${FFTW3_LIBRARIES})
    endif(OPENMP_FOUND)
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
    list(APPEND extra_includes ${FFTW2_INCLUDE_DIR})
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

set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -DBOOST_ENABLE_ASSERTS")

if(ENABLE_WARNINGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
endif(ENABLE_WARNINGS)

if(WARNINGS_ARE_ERRORS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
endif(WARNINGS_ARE_ERRORS)

if(CLANG_FIX)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
endif(CLANG_FIX)

# GSV vector flags
if (NOT VECTOR_FLAGS)
    set(VECTOR_FLAGS "SSE")
endif()

# Additional CXXFLAGS
if(EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif(EXTRA_CXX_FLAGS)

message("-- VECTOR_FLAGS to ${VECTOR_FLAGS}")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DGSV_${VECTOR_FLAGS}")

# simple timer
if(USE_SIMPLE_TIMER)
    message("-- Simple timer profiling enabled")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_SIMPLE_TIMER")
    if(USE_SIMPLE_TIMER_BARRIER)
        message("-- Simple timer MPI_Barrier enabled")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_SIMPLE_TIMER_BARRIER")
    endif(USE_SIMPLE_TIMER_BARRIER)
    if(USE_SIMPLE_TIMER_MEM)
        message("-- Simple timer memory profiling enabled")
        set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -DUSE_SIMPLE_TIMER_MEM")
    endif(USE_SIMPLE_TIMER_MEM)
endif(USE_SIMPLE_TIMER)

# Additional CXXFLAGS
if(EXTRA_CXX_FLAGS)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${EXTRA_CXX_FLAGS}")
endif(EXTRA_CXX_FLAGS)

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

include_directories(BEFORE ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils/vectorclass)
include_directories(BEFORE ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils)
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

function(add_test_executable name)
#  add_executable(${name} EXCLUDE_FROM_ALL ${ARGN})
  add_executable(${name} ${ARGN})
  if(NOT TARGET test-exes)
    add_custom_target(test-exes)
  endif()
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

smart_include_directories("${extra_includes}")

# kokkos
add_subdirectory(src/synergia/utils/kokkos)
include_directories(${Kokkos_INCLUDE_DIRS_RET})

add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(archived-applications)
add_subdirectory(docs)
add_subdirectory(synergia-script-templates)

