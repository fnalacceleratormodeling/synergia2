add_library(synergia_space_charge SHARED rectangular_grid_domain.cc rectangular_grid.cc
    distributed_rectangular_grid.cc deposit.cc space_charge_3d_open_hockney.cc)
target_link_libraries(synergia_space_charge synergia_bunch
    synergia_simulation synergia_distributed_fft)
#target_link_libraries(bunch synergia_parallel_utils synergia_hdf5_utils foundation)
#add_python_extension(space_charge space_charge_wrap.cc)
#target_link_libraries(space_charge synergia_space_charge ${Boost_LIBRARIES})

install(TARGETS synergia_space_charge DESTINATION lib)
install(FILES
    deposit.h
    distributed_rectangular_grid.h
    rectangular_grid_domain.h
    rectangular_grid.h
    space_charge_3d_open_hockney.h
    DESTINATION include/synergia/collective)

add_subdirectory(tests)
