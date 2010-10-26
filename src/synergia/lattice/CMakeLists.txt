add_library(synergia_lattice SHARED lattice_element.cc lattice_element_slice.cc
    lattice.cc chef_lattice.cc chef_utils.cc element_adaptor.cc)
target_link_libraries(synergia_lattice synergia_foundation ${CHEF_LIBS} ${Boost_REGEX_LIBRARY}
    ${Boost_SERIALIZATION_LIBRARY})

add_python_extension(lattice lattice_wrap.cc)
target_link_libraries(lattice synergia_lattice synergia_foundation ${CHEF_LIBS} ${Boost_LIBRARIES})

install(TARGETS synergia_lattice DESTINATION lib)
install(FILES
    chef_elements.h
    chef_lattice.h
    chef_utils.h
    element_adaptor.h
    lattice_element.h
    lattice_element_slice.h
    lattice.h
    DESTINATION include/synergia/lattice)
install(FILES
    __init__.py
    mad8_parser.py
    mad8_reader.py
    DESTINATION lib/synergia/lattice)
install(TARGETS
    lattice
    DESTINATION lib/synergia/lattice)

add_subdirectory(tests)
