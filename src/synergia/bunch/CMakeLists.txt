file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_library(synergia_bunch fixed_t_z_converter.cc bunch.cc
    bunch_train.cc
    diagnostics.cc core_diagnostics.cc diagnostics_basic.cc
    diagnostics_full2.cc diagnostics_particles.cc diagnostics_track.cc
    diagnostics_reference_particle.cc diagnostics_bulk_track.cc diagnostics_phase_space_density.cc
    populate.cc analysis.cc period.cc diagnostics_apertures_loss.cc)
target_link_libraries(synergia_bunch synergia_parallel_utils synergia_distributed_fft
    synergia_hdf5_utils synergia_foundation synergia_serialization)
if (BUILD_PYTHON_BINDINGS)
    add_python_extension(bunch bunch_wrap.cc)
    target_link_libraries(bunch synergia_bunch synergia_distributed_fft ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_bunch DESTINATION lib)
install(FILES
    bunch.h
    bunch_train.h
    core_diagnostics.h
    diagnostics.h
    diagnostics_basic.h
    diagnostics_full2.h
    diagnostics_track.h
    diagnostics_bulk_track.h
    diagnostics_phase_space_density.h
    diagnostics_particles.h
    diagnostics_apertures_loss.h
    fixed_t_z_converter.h
    populate.h
    analysis.h
    period.h
    DESTINATION include/synergia/bunch)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION lib/synergia/bunch)
    install(TARGETS
        bunch
        DESTINATION lib/synergia/bunch)
endif ()

add_subdirectory(tests)
