add_library(bunch SHARED fixed_t_z_converter.cc bunch.cc diagnostics.cc
    diagnostics_writer.cc populate.cc)
target_link_libraries(bunch parallel_utils hdf5_utils foundation)
add_python_extension(pybunch bunch_wrap.cc)
target_link_libraries(pybunch bunch ${Boost_LIBRARIES})

install(TARGETS bunch DESTINATION lib)
install(FILES bunch.h
    diagnostics.h
    diagnostics_writer.h
    fixed_t_z_converter.h
    populate.h
    DESTINATION include/synergia/bunch)

add_subdirectory(tests)
