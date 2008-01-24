project(synergia2)

set(CMAKE_VERBOSE_MAKEFILE ON)

INCLUDE(FindPythonInterp)

INCLUDE(FindPythonLibs)
INCLUDE_DIRECTORIES(${PYTHON_INCLUDE_PATH})
LINK_DIRECTORIES(${PYTHON_LIBRARIES})

FIND_PACKAGE(Boost)
INCLUDE_DIRECTORIES(${Boost_INCLUDE_DIRS})
LINK_DIRECTORIES(${Boost_LIBRARY_DIRS})
set(BOOST_LIBS -lboost_python-gcc)

set(CHEF_INSTALL_DIR /home2/amundson/work/synergia2-base/install)
set(CHEF_LIBS -lphysics_toolkit -lbeamline -lbmlfactory -lglib-2.0 -lbasic_toolkit -lmxyzptlk)
INCLUDE_DIRECTORIES(${CHEF_INSTALL_DIR}/include)
LINK_DIRECTORIES(${CHEF_INSTALL_DIR}/lib)


# add_subdirectory(s2_fish)

add_library(apply_map MODULE apply_map.cc)
target_link_libraries(apply_map  ${CHEF_LIBS} ${BOOST_LIBS})
set_target_properties(apply_map PROPERTIES PREFIX "")

add_library(chef_propagate MODULE chef_propagate.cc)
target_link_libraries(chef_propagate  ${CHEF_LIBS} ${BOOST_LIBS})
set_target_properties(chef_propagate PROPERTIES PREFIX "")

add_library(error_eater MODULE error_eater.cc)
target_link_libraries(error_eater ${BOOST_LIBS})
set_target_properties(error_eater PROPERTIES PREFIX "")

add_library(mappers MODULE mappers.cc)
target_link_libraries(mappers  ${CHEF_LIBS} ${BOOST_LIBS})
set_target_properties(mappers PROPERTIES PREFIX "")

add_library(hist2d MODULE hist2d.cc)
target_link_libraries(hist2d ${BOOST_LIBS})
set_target_properties(hist2d PROPERTIES PREFIX "")
