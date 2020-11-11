#file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
#    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

if(${Kokkos_ENABLE_CUDA})
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
    #diagnostics_basic.cc
    diagnostics_full2.cc 
    diagnostics_full2_host.cc 
    diagnostics_particles.cc 
    #diagnostics_track.cc
    diagnostics_loss.cc
    #diagnostics_reference_particle.cc 
    diagnostics_bulk_track.cc 
    #diagnostics_bulk_spectator_track.cc 
    #diagnostics_phase_space_density.cc
    populate.cc 
    populate_global.cc 
    populate_host.cc 
    #analysis.cc 
    period.cc 
    )

target_link_libraries(synergia_bunch 
    synergia_parallel_utils 
    synergia_hdf5_utils 
    synergia_foundation 
    synergia_serialization 
    kokkos
    )

if (BUILD_PYTHON_BINDINGS)
    #add_python_extension(bunch bunch_wrap.cc)
    #target_link_libraries(bunch synergia_bunch synergia_distributed_fft ${Boost_LIBRARIES})

    add_library(synergia_pydiag diagnostics_py.cc)

    # the -fvisibility flag suppress the warning, but it also
    # prevents the pydiag from exposing to Python as well
    #target_compile_options(synergia_pydiag
    #    PRIVATE "-fvisibility=hidden")

    target_include_directories(synergia_pydiag 
        PRIVATE ${PYBIND11_INCLUDE_DIR}
        PRIVATE ${pybind11_INCLUDE_DIR}
        PRIVATE ${PYTHON_INCLUDE_DIRS} )

    target_link_libraries(synergia_pydiag 
        PRIVATE ${PYTHON_LIBRARIES})

    # add pydiag to synergia_bunch library
    target_link_libraries(synergia_bunch synergia_pydiag)

    pybind11_add_module(bunch bunch_pywrap.cc)
    target_link_libraries(bunch PRIVATE synergia_bunch)

    copy_file(__init__.py bunch)
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
    diagnostics_bulk_spectator_track.h
    diagnostics_phase_space_density.h
    diagnostics_particles.h
    fixed_t_z_converter.h
    populate.h
    analysis.h
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
