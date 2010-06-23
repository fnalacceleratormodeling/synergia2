add_library(space_charge SHARED rectangular_grid_domain.cc rectangular_grid.cc
    distributed_rectangular_grid.cc deposit.cc)
target_link_libraries(space_charge bunch)
#target_link_libraries(bunch parallel_utils hdf5_utils foundation)
#add_python_extension(pyspace_charge space_charge_wrap.cc)
#target_link_libraries(pyspace_charge space_charge ${Boost_LIBRARIES})

add_subdirectory(tests)
