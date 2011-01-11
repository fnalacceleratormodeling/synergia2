add_library(synergia_bunch SHARED fixed_t_z_converter.cc bunch.cc
    diagnostics.cc diagnostics_writer.cc populate.cc)
# Temporary fix for Fedora 14
set_source_files_properties(populate.cc PROPERTIES COMPILE_FLAGS -O0)
# End temporary fix
target_link_libraries(synergia_bunch synergia_parallel_utils
    synergia_hdf5_utils synergia_foundation)
add_python_extension(bunch bunch_wrap.cc)
target_link_libraries(bunch synergia_bunch ${Boost_LIBRARIES})

install(TARGETS synergia_bunch DESTINATION lib)
install(FILES bunch.h
    diagnostics.h
    diagnostics_writer.h
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
