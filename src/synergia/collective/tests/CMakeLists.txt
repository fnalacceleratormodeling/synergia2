add_executable(test_deposit test_deposit.cc)
target_link_libraries(test_deposit synergia_collective synergia_test_main)
add_mpi_test(test_deposit 1)

add_executable(test_space_charge_3d_open_hockney_mpi
               test_space_charge_3d_open_hockney_mpi.cc)
target_link_libraries(test_space_charge_3d_open_hockney_mpi synergia_collective
                      synergia_serialization synergia_test_main)
add_mpi_test(test_space_charge_3d_open_hockney_mpi 1)

add_executable(test_space_charge_3d_rectangular_mpi
               test_space_charge_3d_rectangular_mpi.cc)
target_link_libraries(test_space_charge_3d_rectangular_mpi synergia_collective
                      synergia_serialization synergia_test_main)
add_mpi_test(test_space_charge_3d_rectangular_mpi 1)

if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_impedance_options.py)
endif()

if(${BUILD_FD_SPACE_CHARGE_SOLVER})
  add_executable(test_space_charge_3d_fd_mpi test_space_charge_3d_fd_mpi.cc)
  target_link_libraries(test_space_charge_3d_fd_mpi synergia_collective
                        synergia_serialization synergia_test_main)
  add_mpi_test(test_space_charge_3d_fd_mpi 1)

endif()
