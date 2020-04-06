cmake_minimum_required(VERSION 3.10 FATAL_ERROR)
project(SYNERGIA2 VERSION 2.0.90 LANGUAGES CXX C)

enable_testing()
set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(CMAKE_MODULE_PATH "${SYNERGIA2_SOURCE_DIR}/CMake")

if ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # using Clang
  set(COMPILER_SPECIFIC_CXX_FLAGS "-Wno-#pragma-messages -Wno-potentially-evalauted-expression -Wno-register")
  message("Using clang++")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "AppleClang")
    # using Apples's clang
    set(COMPILER_SPECIFIC_CXX_FLAGS "-Wno-#pragma-messages -Wno-potentially-evaluated-expression -Wno-register")
    message("Using Apple clang++")
elseif ("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  # using GCC
  set(COMPILER_SPECIFIC_CXX_FLAGS "-Wno-error=unused-result -Wno-register")
  set(COMPILER_SPECIFIC_SHARED_LINKER_FLAGS "-Wl,-z,defs")
  message("Using g++")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Intel")
  # using Intel C++
  set(COMPILER_SPECIFIC_CXX_FLAGS "-Wno-error=unused-result")
  set(COMPILER_SPECIFIC_SHARED_LINKER_FLAGS "-Wl,-z,defs")
  message("Using Intel C++")
else()
  message(WARNING "I can not tell what compiler is being used. Will proceed with generic settings.")
endif()

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall ${COMPILER_SPECIFIC_CXX_FLAGS} -Werror")
set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} ${COMPILER_SPECIFIC_SHARED_LINKER_FLAGS}")
set(CMAKE_MODULE_LINKER_FLAGS "${CMAKE_MODULE_LINKER_FLAGS} ${COMPILER_SPECIFIC_SHARED_LINKER_FLAGS}")

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
# needed to avoid warnings about rpath policy on Mac OSX
if(POLICY CMP0042)
    cmake_policy(SET CMP0042 NEW)
endif()

# get version number generation right
if(POLICY CMP0048)
    cmake_policy(SET CMP0048 NEW)
endif()

# use <Package>_ROOT variables
if(POLICY CMP0074)
  cmake_policy(SET CMP0074 NEW)
endif()

##
# Find necessary packages
##

# python
if (USE_PYTHON_3)
  message(STATUS "Trying to find python 3")
  message(STATUS "A failure here is not fatal... we'll try differently later")
  set(Python_ADDITIONAL_VERSIONS "3")
  find_package(Python3 COMPONENTS Interpreter Development)
  if (Python3_FOUND)
    message(STATUS "We found Python 3")
    message(STATUS "Python3_LIBRARIES: ${Python3_LIBRARIES}")
    set(MY_PYTHON_LIBRARY ${Python3_LIBRARIES})
    set(MY_PYTHON_EXECUTABLE ${Python3_EXECUTABLE})
    set(MY_PYTHON_INCLUDE_DIRECTORY ${Python3_INCLUDE_DIRS})
    set(MY_PYTHON_VERSION_MAJOR ${Python3_VERSION_MAJOR})
    set(MY_PYTHON_VERSION_MINOR ${Python3_VERSION_MINOR})
    # If we could use cmake v3.12 or newer, this is simple:
    #find_package(Python3 COMPONENTS Interpreter Development)
  endif()
else()
  message(STATUS "Trying to find python 2")
  message(STATUS "A failure here is not fatal... we'll try differently later")
  set(Python_ADDTIONAL_VERSIONS "2")
  find_package(Python2 COMPONENTS Interpreter Development)
  if (Python2_FOUND)
    set(MY_PYTHON_LIBRARY ${Python2_LIBRARIES})
    set(MY_PYTHON_EXECUTABLE ${Python2_EXECUTABLE})
    set(MY_PYTHON_INCLUDE_DIRECTORY ${Python2_INCLUDE_DIRS})
    set(MY_PYTHON_VERSION_MAJOR ${Python2_VERSION_MAJOR})
    set(MY_PYTHON_VERSION_MINOR ${Python2_VERSION_MINOR})
    # If we could use cmake v3.12 or newer, this is simple:
    #find_package(Python3 COMPONENTS Interpreter Development)
  endif()
endif()

# If we could use cmake v3.12 or newer, this would not be needed.
if (NOT MY_PYTHON_LIBRARY)
  find_package(PythonInterp REQUIRED)
  find_package(PythonLibs REQUIRED)
  set(MY_PYTHON_LIBRARY ${PYTHON_LIBRARIES})
  set(MY_PYTHON_EXECUTABLE ${PYTHON_EXECUTABLE})
  set(MY_PYTHON_INCLUDE_DIRECTORY ${PYTHON_INCLUDE_DIRS})
  set(MY_PYTHON_VERSION_MAJOR ${PYTHON_VERSION_MAJOR})
  set(MY_PYTHON_VERSION_MINOR ${PYTHON_VERSION_MINOR})
endif()
include_directories(${MY_PYTHON_INCLUDE_DIR})

message(STATUS "Done looking for Python")
message(STATUS "MY_PYTHON_VERSION_MAJOR is: ${MY_PYTHON_VERSION_MAJOR}")
message(STATUS "MY_PYTHON_VERSION_MINOR is: ${MY_PYTHON_VERSION_MINOR}")

include(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)

# boost
set(Boost_NO_BOOST_CMAKE ON) # Do *not* use CMake support from Boost.
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREAD ON)
find_package(Boost
             REQUIRED
             COMPONENTS regex unit_test_framework serialization system filesystem)
set(first_boost_libraries ${Boost_LIBRARIES})
set(first_boost_library_dirs ${Boost_LIBRARY_DIRS})

if (USE_PYTHON_3)
    set(pstem python3 python35 python36 python37)
else()
    set(pstem python python27)
endif()

# Go through all relevant possible choices of boost library name, stopping
# when we find one that works.
foreach (plib ${pstem})
    find_package(Boost QUIET COMPONENTS ${plib})
    if (Boost_LIBRARIES)
      string(TOUPPER "${plib}" MY_PYLIB)
      set(MY_BOOST_PYTHON_LIBRARY "${Boost_${MY_PYLIB}_LIBRARY}")
      break()
    endif()
endforeach()
set(Boost_LIBRARIES ${first_boost_libraries} ${Boost_LIBRARIES})
set(Boost_LIBRARY_DIRS ${first_boost_library_dirs} ${Boost_LIBRARY_DIRS})

include_directories(${Boost_INCLUDE_DIRS})
link_directories(${Boost_LIBRARY_DIRS})

message(STATUS "Boost_INCLUDE_DIRS = ${Boost_INCLUDE_DIRS}")
message(STATUS "Boost_LIBRARY_DIRS = ${Boost_LIBRARY_DIRS}")
message(STATUS "Boost_LIBRARIES    = ${Boost_LIBRARIES}")
message(STATUS "Boost_SERIALIZATION_LIBRARY = ${Boost_SERIALIZATION_LIBRARY}")
message(STATUS "MY_BOOST_PYTHON_LIBRARY = ${MY_BOOST_PYTHON_LIBRARY}")

# MPI
find_package(MPI REQUIRED CXX)
execute_process(COMMAND ${MPIEXEC_EXECUTABLE} --version OUTPUT_VARIABLE MPI_VERSION)
message(STATUS "MPIEXEC_EXECUTABLE is ${MPIEXEC_EXECUTABLE}")
if(${MPI_VERSION} MATCHES "OpenRTE.* ")
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

# hdf5
find_package(HDF5 REQUIRED)
include_directories(${HDF5_INCLUDE_DIRS})
message(STATUS "${HDF5_LIBRARIES}")

# gsl
find_package(GSL REQUIRED)
include_directories(${GSL_INCLUDE_DIRS})
message(STATUS "Found GSL at ${GSL_INCLUDE_DIRS}")

# chef
find_package(CHEF REQUIRED)
include_directories(${CHEF_INCLUDE_DIR})
link_directories(${CHEF_LIB_DIR})

# numpy
find_package(NUMPY REQUIRED)
include_directories(${NUMPY_INCLUDE_DIR})

# mpi4py
find_package(MPI4PY REQUIRED)
include_directories(${MPI4PY_INCLUDE_DIR})

# misc
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)

# Eigen3
find_package(Eigen3 REQUIRED)
include_directories(${EIGEN3_INCLUDE_DIR})
#list(APPEND extra_includes ${EIGEN3_INCLUDE_DIR})

# fftw3 (and, possibly, fftw2)
find_package(FFTW3 REQUIRED)
message(STATUS "FFTW3 was found")
message(STATUS "FFTW3_FOUND: ${FFTW3_FOUND}")
message(STATUS "FFTW3_MPI_LIBRARIES: ${FFTW3_MPI_LIBRARIES}")
if(FFTW3_MPI_FOUND)
    include_directories(${FFTW3_INCLUDE_DIR})
    #list(APPEND extra_includes ${FFTW3_INCLUDE_DIR})
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
    include_directories(${FFTW2_INCLUDE_DIR})
    #list(APPEND extra_includes ${FFTW2_INCLUDE_DIR})
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
else()
  message(STATUS "Simple time profiling not enabled")
endif()

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
if ( (NOT BUILD_SHARED_LIBS))
  message(WARNING "Static library build selected, ignored for python bindings")
endif()

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
set(PYTHON_INSTALL_DIR "lib/python${MY_PYTHON_VERSION_MAJOR}.${MY_PYTHON_VERSION_MINOR}/site-packages/"
    CACHE PATH "python install directory")

set(SYNERGIA_EXTRA_LIBRARY_PATHS "" CACHE STRING "additional library paths for synergia executable")
set(SYNERGIA_EXTRA_PYTHON_PATHS "" CACHE STRING "additional python module paths for synergia executable")

set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH "${CMAKE_INSTALL_PREFIX}/${LIB_INSTALL_DIR}")
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)

#smart_include_directories("${extra_includes}")
add_subdirectory(src)
add_subdirectory(examples)
add_subdirectory(archived-applications)
add_subdirectory(docs)
add_subdirectory(synergia-script-templates)

