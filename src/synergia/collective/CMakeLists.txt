add_library(synergia_collective rectangular_grid_domain.cc rectangular_grid.cc
    distributed_rectangular_grid.cc fftw_rectangular_helper.cc deposit.cc space_charge_3d_open_hockney.cc
    space_charge_2d_open_hockney.cc space_charge_2d_bassetti_erskine.cc
    space_charge_rectangular.cc impedance.cc wake_field.cc)
target_link_libraries(synergia_collective synergia_bunch
    synergia_simulation synergia_distributed_fft)

if (BUILD_PYTHON_BINDINGS)
    add_python_extension(collective collective_wrap.cc)
    target_link_libraries(collective synergia_collective synergia_bunch
        synergia_simulation synergia_distributed_fft ${Boost_LIBRARIES})
endif ()

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/__init__.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

install(TARGETS synergia_collective DESTINATION lib)
install(FILES
    deposit.h
    distributed_rectangular_grid.h
    interpolate_rectangular_zyx.h
    rectangular_grid_domain.h
    rectangular_grid.h
    fftw_rectangular_helper.h
    space_charge_3d_open_hockney.h
    space_charge_2d_open_hockney.h
    space_charge_rectangular.h
    impedance.h
    wake_field.h
    DESTINATION include/synergia/collective)
if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION lib/synergia/collective)
    install(TARGETS
        collective
        DESTINATION lib/synergia/collective)
endif ()

add_subdirectory(tests)
