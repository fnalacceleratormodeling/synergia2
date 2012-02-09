add_library(synergia_bunch fixed_t_z_converter.cc bunch.cc
    diagnostics.cc multi_diagnostics.cc  diagnostics_write_helper.cc 
        bunch_with_diagnostics.cc train.cc populate.c
	analysis.cc period.cc)
target_link_libraries(synergia_bunch synergia_parallel_utils synergia_distributed_fft
     synergia_hdf5_utils synergia_foundation)
if (BUILD_PYTHON_BINDINGS)
    add_python_extension(bunch bunch_wrap.cc)
    target_link_libraries(bunch synergia_bunch synergia_distributed_fft ${Boost_LIBRARIES})
endif ()

install(TARGETS synergia_bunch DESTINATION lib)
install(FILES bunch.h
    diagnostics.h
    bunch_with_diagnostics.h
    train.h
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
