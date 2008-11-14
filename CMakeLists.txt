project(SYNERGIA2)
cmake_minimum_required(VERSION 2.6)
enable_testing()

INCLUDE(FindPythonInterp)

INCLUDE(FindPythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
LINK_DIRECTORIES(${PYTHON_LIBRARIES})


SET(Boost_USE_STATIC_LIBS OFF)
SET(Boost_USE_MULTITHREAD ON)
FIND_PACKAGE(Boost COMPONENTS python regex unit_test_framework)

IF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 2.6)
    set(Boost_LIBRARIES -lboost_python-gcc)
ENDIF (${CMAKE_MAJOR_VERSION}.${CMAKE_MINOR_VERSION} LESS 2.6)

message("Boost_FOUND = '${Boost_FOUND}'")
message("Boost_ERROR_REASON = '${Boost_ERROR_REASON}'")
message("Boost_INCLUDE_DIRS = '${Boost_INCLUDE_DIRS}'")
message("Boost_INCLUDE_DIR = '${Boost_INCLUDE_DIR}'")
message("Boost_LIBRARIES = '${Boost_LIBRARIES}'")
message("Boost_LIBRARY_DIRS = '${Boost_LIBRARY_DIRS}'")
message("Boost_VERSION = '${Boost_VERSION}'")
message("Boost_LIB_VERSION = '${Boost_LIB_VERSION}'")
message("Boost_MAJOR_VERSION = '${Boost_MAJOR_VERSION}'")
message("Boost_MINOR_VERSION = '${Boost_MINOR_VERSION}'")
message("Boost_SUBMINOR_VERSION = '${Boost_SUBMINOR_VERSION}'")

INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})

INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindCHEF.cmake)
INCLUDE_DIRECTORIES(${CHEF_INCLUDE_DIRS})
LINK_DIRECTORIES(${CHEF_LIBRARY_DIRS})

INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/AddPythonExtension.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindGSL.cmake)
INCLUDE(${SYNERGIA2_SOURCE_DIR}/CMake/FindFFTW2.cmake)

configure_file("${SYNERGIA2_SOURCE_DIR}/synergia-local.in" 
    "${SYNERGIA2_BINARY_DIR}/synergia-local" IMMEDIATE)

include_directories(${SYNERGIA2_SOURCE_DIR})

add_subdirectory(array_nd)
add_subdirectory(synergia)
add_subdirectory(s2_fish)
add_subdirectory(impact)

add_subdirectory(octapy)
