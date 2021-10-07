if(${Kokkos_ENABLE_CUDA})
    set_source_files_properties(mx_parse.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(mx_tree.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(madx.cc PROPERTIES COMPILE_FLAGS --host-only)
endif()

add_library(synergia_lattice 
    lattice_element.cc 
    lattice_element_processor.cc 
    lattice.cc
    lattice_element_slice.cc
    madx.cc 
    mx_expr.cc  
    mx_parse.cc  
    mx_tree.cc
    madx_reader.cc 
    )

target_link_libraries(synergia_lattice 
    synergia_foundation 
    synergia_bunch
    )

target_link_options(synergia_lattice
    PRIVATE ${LINKER_OPTIONS}
    )

if (BUILD_PYTHON_BINDINGS)

    pybind11_add_module(lattice MODULE NO_EXTRAS lattice_pywrap.cc)
    target_link_libraries(lattice PRIVATE synergia_lattice)

    copy_file(__init__.py lattice)
    copy_file(mad8_parser.py lattice)
    copy_file(mad8_reader.py lattice)
    copy_file(simplify.py lattice)

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

    #install(
    #    TARGETS lattice
    #    DESTINATION ${PYTHON_INSTALL_DIR}/synergia/lattice )

endif ()

add_subdirectory(tests)
