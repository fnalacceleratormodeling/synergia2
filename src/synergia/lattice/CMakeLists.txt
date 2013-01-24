add_library(synergia_lattice lattice_element.cc lattice_element_slice.cc
    lattice.cc chef_lattice.cc chef_lattice_section.cc chef_utils.cc
    element_adaptor.cc mad8_adaptors.cc madx_adaptors.cc
    element_adaptor_map.cc mad8_adaptor_map.cc madx_adaptor_map.cc
    lattice_diagnostics.cc
    madx.cc mx_expr.cc  mx_parse.cc  mx_tree.cc
    madx_reader.cc)

# jfa: temporary entry
add_executable(thetest test.cc)
target_link_libraries(thetest synergia_foundation synergia_lattice)

target_link_libraries(synergia_lattice synergia_foundation
    ${CHEF_LIBS} ${Boost_REGEX_LIBRARY} ${Boost_SERIALIZATION_LIBRARY})

if (BUILD_PYTHON_BINDINGS)
    add_python_extension(lattice lattice_wrap.cc)
    target_link_libraries(lattice synergia_lattice synergia_foundation
        ${CHEF_LIBS} ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_lattice DESTINATION lib)
install(FILES
    chef_elements.h
    chef_lattice.h
    chef_lattice_section.h
    chef_lattice_section_fwd.h
    chef_utils.h
    element_adaptor.h
    lattice_element.h
    lattice_element_slice.h
    lattice.h
    madx.h
    madx_reader.h
    mx_expr.h
    mx_parse.h
    mx_tree.h
    DESTINATION include/synergia/lattice)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        mad8_parser.py
        mad8_reader.py
        simplify.py
        DESTINATION lib/synergia/lattice)
    install(TARGETS
        lattice
        DESTINATION lib/synergia/lattice)
endif ()

add_subdirectory(tests)
