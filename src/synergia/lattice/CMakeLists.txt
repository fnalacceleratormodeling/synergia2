file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mad8_parser.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/mad8_reader.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/simplify.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

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
    #${CHEF_LIBS} 
    #${Boost_REGEX_LIBRARY} 
    #${Boost_SERIALIZATION_LIBRARY}
    )

if (BUILD_PYTHON_BINDINGS)
    #add_python_extension(lattice lattice_wrap.cc)
    #target_link_libraries(lattice synergia_lattice synergia_foundation synergia_bunch
    #    ${CHEF_LIBS} ${Boost_LIBRARIES})

    pybind11_add_module(lattice MODULE NO_EXTRAS lattice_pywrap.cc)
    target_link_libraries(lattice PRIVATE synergia_lattice)

    copy_file(__init__.py lattice)
    copy_file(mad8_parser.py lattice)
    copy_file(mad8_reader.py lattice)

endif ()

install(TARGETS synergia_lattice DESTINATION ${LIB_INSTALL_DIR})

install(FILES
    lattice_element.h
    lattice_element_slice.h
    lattice.h
    madx.h
    madx_reader.h
    mx_expr.h
    mx_parse.h
    mx_tree.h
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
