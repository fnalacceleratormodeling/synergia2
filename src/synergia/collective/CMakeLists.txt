add_library(synergia_collective 
    deposit.cc 
    space_charge_3d_open_hockney.cc
    space_charge_2d_open_hockney.cc 
    space_charge_2d_kv.cc
    space_charge_rectangular.cc 
    impedance.cc 
    wake_field.cc 
    )

target_link_libraries(synergia_collective 
    synergia_bunch
    synergia_distributed_fft
    synergia_parallel_utils
    synergia_bunchsim
    synergia_serialization
    )

target_link_options(synergia_collective
    PRIVATE ${LINKER_OPTIONS}
    )

if (BUILD_PYTHON_BINDINGS)
    pybind11_add_module(collective MODULE NO_EXTRAS collective_pywrap.cc)
    target_link_libraries(collective PRIVATE synergia_collective)
    copy_file(__init__.py collective)
endif ()

install(TARGETS synergia_collective DESTINATION lib)
install(FILES
    deposit.h
    rectangular_grid_domain.h
    rectangular_grid.h
    space_charge_3d_open_hockney.h
    space_charge_2d_open_hockney.h
    space_charge_2d_kv.h
    space_charge_rectangular.h
    impedance.h
    wake_field.h
    DESTINATION ${INCLUDE_INSTALL_DIR}/synergia/collective)

if (BUILD_PYTHON_BINDINGS)
    install(FILES
        __init__.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/collective)
    install(TARGETS
        collective
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia/collective)

endif ()

add_subdirectory(tests)
