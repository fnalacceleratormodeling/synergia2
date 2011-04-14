add_library(synergia_bunch SHARED fixed_t_z_converter.cc bunch.cc
    diagnostics.cc multi_diagnostics.cc diagnostics_write_helper.cc populate.cc)
target_link_libraries(synergia_bunch synergia_parallel_utils
    synergia_hdf5_utils synergia_foundation)
add_python_extension(bunch bunch_wrap.cc)
target_link_libraries(bunch synergia_bunch ${Boost_LIBRARIES})

install(TARGETS synergia_bunch DESTINATION lib)
install(FILES bunch.h
    diagnostics.h
    multi_diagnostics.h
    fixed_t_z_converter.h
    populate.h
    DESTINATION include/synergia/bunch)
install(FILES
    __init__.py
    DESTINATION lib/synergia/bunch)
install(TARGETS
    bunch
    DESTINATION lib/synergia/bunch)

add_subdirectory(tests)
