if((NOT ${USE_EXTERNAL_KOKKOS}) AND (DEFINED ${Kokkos_ENABLE_CUDA}))
    # We are using nvcc-wrapper, need to set these flags for now.	
    set_source_files_properties(populate_host.cc PROPERTIES COMPILE_FLAGS --host-only)
    set_source_files_properties(diagnostics_full2_host.cc PROPERTIES COMPILE_FLAGS --host-only)
endif()

add_library(synergia_bunch 
    bunch_particles.cc
    bunch.cc
    bunch_train.cc
    core_diagnostics.cc 
    diagnostics_worker.cc 
    diagnostics_file.cc 
    diagnostics_full2.cc 
    diagnostics_full2_host.cc 
    diagnostics_particles.cc 
    diagnostics_loss.cc
    diagnostics_bulk_track.cc 
    populate.cc 
    populate_global.cc 
    populate_host.cc 
    period.cc 
    )

target_link_libraries(synergia_bunch 
    synergia_parallel_utils 
    synergia_hdf5_utils 
    synergia_foundation 
    synergia_serialization 
    ${kokkos_libs}
    )

target_link_options(synergia_bunch
    PRIVATE ${LINKER_OPTIONS}
    )

if (BUILD_PYTHON_BINDINGS)
    add_library(synergia_pydiag diagnostics_py.cc)

    # the -fvisibility flag suppress the warning, but it also
    # prevents the pydiag from exposing to Python as well
    #set_target_properties(synergia_pydiag 
    #    PROPERTIES CXX_VISIBILITY_PRESET "hidden")

    target_link_libraries(synergia_pydiag
        PRIVATE synergia_bunch pybind11::module)

    pybind11_add_module(bunch MODULE NO_EXTRAS bunch_pywrap.cc)
    target_link_libraries(bunch PRIVATE synergia_bunch synergia_pydiag)

    copy_file(__init__.py bunch)

    install(FILES
            __init__.py
            DESTINATION ${PYTHON_INSTALL_DIR}/synergia/bunch)
    install(TARGETS
            bunch
            synergia_pydiag
            DESTINATION ${LIB_INSTALL_DIR})

endif ()


install(TARGETS synergia_bunch DESTINATION ${LIB_INSTALL_DIR})
install(FILES
    bunch.h
    bunch_train.h
    core_diagnostics.h
    diagnostics.h
    diagnostics_basic.h
    diagnostics_full2.h
    diagnostics_track.h
    diagnostics_bulk_track.h
    diagnostics_particles.h
    fixed_t_z_converter.h
    populate.h
    period.h
    DESTINATION ${INCLUDE_INSTALL_DIR}/synergia/bunch)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/bunch)
    install(TARGETS
        bunch
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/bunch)
endif ()

add_subdirectory(tests)
