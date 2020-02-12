
add_executable(test_hdf5_file_mpi test_hdf5_file_mpi.cc)
target_link_libraries(test_hdf5_file_mpi synergia_hdf5_utils synergia_test_main)

add_mpi_test(test_hdf5_file_mpi 1)
add_mpi_test(test_hdf5_file_mpi 2)
add_mpi_test(test_hdf5_file_mpi 3)
add_mpi_test(test_hdf5_file_mpi 4)

#add_executable(test_bunch_simulator test_bunch_simulator.cc)
#target_link_libraries(test_bunch_simulator synergia_simulation synergia_test_main)
#add_test(NAME test_bunch_simulator COMMAND test_bunch_simulator)



