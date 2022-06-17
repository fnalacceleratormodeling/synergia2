add_executable(
    test_deposit
    test_deposit.cc)
target_link_libraries(
    test_deposit 
    synergia_collective
    synergia_test_main
    )
add_mpi_test(test_deposit 1)


add_executable(
    test_space_charge_3d_open_hockney_mpi 
    test_space_charge_3d_open_hockney_mpi.cc)
target_link_libraries(
    test_space_charge_3d_open_hockney_mpi 
    synergia_collective
    synergia_serialization
    synergia_test_main
    )
add_mpi_test(test_space_charge_3d_open_hockney_mpi 1)
