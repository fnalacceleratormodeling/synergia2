add_library(synergia_collective SHARED rectangular_grid_domain.cc rectangular_grid.cc
    distributed_rectangular_grid.cc deposit.cc space_charge_3d_open_hockney.cc
    space_charge_2d_bassetti_erskine.cc)
target_link_libraries(synergia_collective synergia_bunch
    synergia_simulation synergia_distributed_fft)

add_python_extension(collective collective_wrap.cc)
target_link_libraries(collective synergia_collective synergia_bunch
    synergia_simulation synergia_distributed_fft ${Boost_LIBRARIES})


install(TARGETS synergia_collective DESTINATION lib)
install(FILES
    deposit.h
    distributed_rectangular_grid.h
    rectangular_grid_domain.h
    rectangular_grid.h
    space_charge_3d_open_hockney.h
    DESTINATION include/synergia/collective)
install(FILES
    __init__.py
    DESTINATION lib/synergia/collective)
install(TARGETS
    collective
    DESTINATION lib/synergia/collective)

add_subdirectory(tests)
