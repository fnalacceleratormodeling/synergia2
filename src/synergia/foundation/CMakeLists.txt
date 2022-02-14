
add_library(synergia_foundation 
    four_momentum.cc 
    reference_particle.cc
    distribution.cc 
    )

target_link_libraries(synergia_foundation 
    ${GSL_CONFIG_LIBS}
    synergia_parallel_utils 
    lsexpr
    gsl
    gslcblas
    ${kokkos_libs}
    )

target_link_options(synergia_foundation
    PRIVATE ${LINKER_OPTIONS}
    )

if (BUILD_PYTHON_BINDINGS)

    if(${Kokkos_ENABLE_CUDA})
        #set_source_files_properties(foundation_pywrap.cc 
        #    PROPERTIES COMPILE_FLAGS --host-only)
    endif()

    pybind11_add_module(foundation MODULE NO_EXTRAS foundation_pywrap.cc)
    target_link_libraries(foundation PRIVATE synergia_foundation)

    copy_file(__init__.py foundation)
endif ()

install(TARGETS synergia_foundation DESTINATION ${LIB_INSTALL_DIR})
install(FILES
    distribution.h
    #diagnostics_write_helper.h
    four_momentum.h
    math_constants.h
    #multi_diagnostics.h
    physical_constants.h
    reference_particle.h
    DESTINATION ${INCLUDE_INSTALL_DIR}/synergia/foundation)

if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/foundation)
    install(TARGETS
        foundation
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/foundation)
endif ()

add_subdirectory(tests)
