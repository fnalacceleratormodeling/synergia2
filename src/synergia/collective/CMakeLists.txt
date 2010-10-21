add_library(space_charge SHARED rectangular_grid_domain.cc rectangular_grid.cc
    distributed_rectangular_grid.cc deposit.cc space_charge_3d_open_hockney.cc)
target_link_libraries(space_charge bunch simulation distributed_fft)
#target_link_libraries(bunch parallel_utils hdf5_utils foundation)
#add_python_extension(pyspace_charge space_charge_wrap.cc)
#target_link_libraries(pyspace_charge space_charge ${Boost_LIBRARIES})

install(TARGETS space_charge DESTINATION lib)
install(FILES
    deposit.h
    distributed_rectangular_grid.h
    rectangular_grid_domain.h
    rectangular_grid.h
    space_charge_3d_open_hockney.h
    DESTINATION include/synergia/collective)

add_subdirectory(tests)
