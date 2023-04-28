if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_parallel_utils.py)
endif()

add_executable(test_command_line_arg test_command_line_arg.cc)
target_link_libraries(test_command_line_arg synergia_command_line
                      synergia_test_main)
add_mpi_test(test_command_line_arg 1)

if(NOT ${USE_OPENPMD_IO})
  add_executable(test_hdf5_file_mpi test_hdf5_file_mpi.cc)
  target_link_libraries(test_hdf5_file_mpi synergia_hdf5_utils
                        synergia_test_main)

  add_mpi_test(test_hdf5_file_mpi 1)
  add_mpi_test(test_hdf5_file_mpi 2)
  add_mpi_test(test_hdf5_file_mpi 3)
  add_mpi_test(test_hdf5_file_mpi 4)

  add_executable(test_hdf5_write_mpi test_hdf5_write_mpi.cc)
  target_link_libraries(test_hdf5_write_mpi synergia_hdf5_utils
                        synergia_test_main)

  add_mpi_test(test_hdf5_write_mpi 1)
  add_mpi_test(test_hdf5_write_mpi 2)
  add_mpi_test(test_hdf5_write_mpi 3)
  add_mpi_test(test_hdf5_write_mpi 4)

  add_executable(test_hdf5_append_mpi test_hdf5_append_mpi.cc)
  target_link_libraries(test_hdf5_append_mpi synergia_hdf5_utils
                        synergia_test_main)

  add_mpi_test(test_hdf5_append_mpi 1)
  add_mpi_test(test_hdf5_append_mpi 2)
  add_mpi_test(test_hdf5_append_mpi 3)
  add_mpi_test(test_hdf5_append_mpi 4)

  add_executable(test_hdf5_read_mpi test_hdf5_read_mpi.cc)
  target_link_libraries(test_hdf5_read_mpi synergia_hdf5_utils
                        synergia_test_main)

  add_mpi_test(test_hdf5_read_mpi 1)
  add_mpi_test(test_hdf5_read_mpi 2)
  add_mpi_test(test_hdf5_read_mpi 3)
  add_mpi_test(test_hdf5_read_mpi 4)
endif()

if("${ENABLE_KOKKOS_BACKEND}" STREQUAL "CUDA")
  add_executable(test_kokkos test_kokkos.cc)
  target_link_libraries(test_kokkos synergia_test_main)

  add_mpi_test(test_kokkos 1)
endif()

add_executable(test_commxx_mpi test_commxx_mpi.cc)
target_link_libraries(test_commxx_mpi synergia_parallel_utils
                      synergia_test_main)
add_mpi_test(test_commxx_mpi 1)
add_mpi_test(test_commxx_mpi 2)
add_mpi_test(test_commxx_mpi 3)
add_mpi_test(test_commxx_mpi 4)

add_executable(test_commxx_serdes test_commxx_serdes.cc)
target_link_libraries(test_commxx_serdes synergia_parallel_utils
                      synergia_test_main)
add_mpi_test(test_commxx_serdes 1)
add_mpi_test(test_commxx_serdes 2)
add_mpi_test(test_commxx_serdes 3)
add_mpi_test(test_commxx_serdes 4)

add_executable(test_distributed_fft2d test_distributed_fft2d.cc)
target_link_libraries(test_distributed_fft2d synergia_distributed_fft
                      synergia_test_main)
add_mpi_test(test_distributed_fft2d 1)
add_mpi_test(test_distributed_fft2d 2)
add_mpi_test(test_distributed_fft2d 3)
add_mpi_test(test_distributed_fft2d 4)

add_executable(test_distributed_fft3d test_distributed_fft3d.cc)
target_link_libraries(test_distributed_fft3d synergia_distributed_fft
                      synergia_test_main)
add_mpi_test(test_distributed_fft3d 1)
add_mpi_test(test_distributed_fft3d 2)
add_mpi_test(test_distributed_fft3d 3)
add_mpi_test(test_distributed_fft3d 4)

if(0)
  add_executable(test_loading_kokkos test_loading_kokkos.cc)
  target_link_libraries(test_loading_kokkos PRIVATE dl)
  target_include_directories(
    test_loading_kokkos
    PRIVATE ${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos
            ${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos/core/src
            ${SYNERGIA2_BINARY_DIR}/src/synergia/utils/kokkos/core/src/desul
            ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils/kokkos/tpls/desul/include
            ${SYNERGIA2_SOURCE_DIR}/src/synergia/utils/kokkos/core/src)
  add_test(test_loading_kokkos test_loading_kokkos)
endif()
