file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mad8_parser.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mad8_reader.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/simplify.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(synergia_lattice lattice_element.cc lattice_element_slice.cc
    lattice.cc chef_lattice.cc chef_lattice_section.cc chef_utils.cc
    element_adaptor.cc mad8_adaptors.cc madx_adaptors.cc
    element_adaptor_map.cc mad8_adaptor_map.cc madx_adaptor_map.cc
    lattice_diagnostics.cc
    madx.cc mx_expr.cc  mx_parse.cc  mx_tree.cc
    madx_reader.cc diagnostics_apertures_loss.cc)

target_link_libraries(synergia_lattice synergia_foundation synergia_bunch
    ${CHEF_LIBS} ${Boost_REGEX_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})

if (BUILD_PYTHON_BINDINGS)
    add_python_extension(lattice lattice_wrap.cc)
    target_link_libraries(lattice synergia_lattice synergia_foundation synergia_bunch
        ${CHEF_LIBS} ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_lattice DESTINATION ${LIB_INSTALL_DIR})
install(FILES
    chef_elements.h
    chef_lattice.h
    chef_lattice_section.h
    chef_lattice_section_fwd.h
    chef_utils.h
    element_adaptor.h
    element_adaptor_map.h
    lattice_element.h
    lattice_element_slice.h
    lattice.h
    mad8_adaptor_map.h
    mad8_adaptors.h
    madx.h
    madx_adaptor_map.h
    madx_adaptors.h
    madx_reader.h
    mx_expr.h
    mx_parse.h
    mx_tree.h
    diagnostics_apertures_loss.h
    DESTINATION ${INCLUDE_INSTALL_DIR}/synergia/lattice)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        mad8_parser.py
        mad8_reader.py
        simplify.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/lattice)
    install(TARGETS
        lattice
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/lattice)
endif ()

add_subdirectory(tests)
