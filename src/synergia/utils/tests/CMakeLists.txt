
add_executable(test_hdf5_file_mpi test_hdf5_file_mpi.cc)
target_link_libraries(test_hdf5_file_mpi synergia_hdf5_utils synergia_test_main)

add_mpi_test(test_hdf5_file_mpi 1)
add_mpi_test(test_hdf5_file_mpi 2)
add_mpi_test(test_hdf5_file_mpi 3)
add_mpi_test(test_hdf5_file_mpi 4)

add_executable(test_hdf5_write_mpi test_hdf5_write_mpi.cc)
target_link_libraries(test_hdf5_write_mpi synergia_hdf5_utils synergia_test_main)

add_mpi_test(test_hdf5_write_mpi 1)
add_mpi_test(test_hdf5_write_mpi 2)
add_mpi_test(test_hdf5_write_mpi 3)
add_mpi_test(test_hdf5_write_mpi 4)

add_executable(test_hdf5_append_mpi test_hdf5_append_mpi.cc)
target_link_libraries(test_hdf5_append_mpi synergia_hdf5_utils synergia_test_main)

add_mpi_test(test_hdf5_append_mpi 1)
add_mpi_test(test_hdf5_append_mpi 2)
add_mpi_test(test_hdf5_append_mpi 3)
add_mpi_test(test_hdf5_append_mpi 4)

add_executable(test_hdf5_read_mpi test_hdf5_read_mpi.cc)
target_link_libraries(test_hdf5_read_mpi synergia_hdf5_utils synergia_test_main)

add_mpi_test(test_hdf5_read_mpi 1)
add_mpi_test(test_hdf5_read_mpi 2)
add_mpi_test(test_hdf5_read_mpi 3)
add_mpi_test(test_hdf5_read_mpi 4)


add_executable(test_kokkos test_kokkos.cc)
target_link_libraries(test_kokkos synergia_test_main)

add_mpi_test(test_kokkos 1)

add_executable(test_commxx_mpi test_commxx_mpi.cc)
target_link_libraries(test_commxx_mpi synergia_parallel_utils synergia_test_main)
add_mpi_test(test_commxx_mpi 1)
add_mpi_test(test_commxx_mpi 2)
add_mpi_test(test_commxx_mpi 3)
add_mpi_test(test_commxx_mpi 4)
