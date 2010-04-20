add_library(lattice SHARED lattice_element.cc lattice_element_slice.cc
    lattice.cc chef_lattice.cc chef_utils.cc element_adaptor.cc)
target_link_libraries(lattice foundation ${CHEF_LIBS} ${Boost_REGEX_LIBRARY})

add_python_extension(pylattice lattice_wrap.cc)
target_link_libraries(pylattice lattice foundation ${CHEF_LIBS} ${Boost_LIBRARIES})

add_subdirectory(tests)
