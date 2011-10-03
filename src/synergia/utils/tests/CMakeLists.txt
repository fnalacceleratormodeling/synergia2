add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_parallel_utils test_parallel_utils.cc)
target_link_libraries(test_parallel_utils synergia_parallel_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_parallel_utils ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_parallel_utils)
add_test(test_parallel_utils_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_parallel_utils.py)

add_executable(test_floating_point test_floating_point.cc)
target_link_libraries(test_floating_point ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_floating_point test_floating_point)

add_executable(test_fast_int_floor test_fast_int_floor.cc)
target_link_libraries(test_fast_int_floor ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_fast_int_floor test_fast_int_floor)

add_executable(test_hdf5_writer test_hdf5_writer.cc)
target_link_libraries(test_hdf5_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_writer test_hdf5_writer)

add_executable(test_hdf5_file test_hdf5_file.cc)
target_link_libraries(test_hdf5_file synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_file test_hdf5_file)

add_executable(test_hdf5_serial_writer test_hdf5_serial_writer.cc)
target_link_libraries(test_hdf5_serial_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_serial_writer test_hdf5_serial_writer)

add_executable(test_hdf5_chunked_array2d_writer test_hdf5_chunked_array2d_writer.cc)
target_link_libraries(test_hdf5_chunked_array2d_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_chunked_array2d_writer test_hdf5_chunked_array2d_writer)

add_executable(test_distributed_fft3d test_distributed_fft3d.cc)
target_link_libraries(test_distributed_fft3d synergia_distributed_fft ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft3d ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_distributed_fft3d)

add_executable(test_distributed_fft3d_mpi test_distributed_fft3d_mpi.cc)
target_link_libraries(test_distributed_fft3d_mpi synergia_distributed_fft
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft3d_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_distributed_fft3d_mpi)
add_test(test_distributed_fft3d_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_distributed_fft3d_mpi)
add_test(test_distributed_fft3d_mpi3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 test_distributed_fft3d_mpi)
add_test(test_distributed_fft3d_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_distributed_fft3d_mpi)

add_executable(test_distributed_fft2d test_distributed_fft2d.cc)
target_link_libraries(test_distributed_fft2d synergia_distributed_fft ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft2d ${MPIEXEC} -np 1 test_distributed_fft2d)

add_executable(test_distributed_fft2d_mpi test_distributed_fft2d_mpi.cc)
target_link_libraries(test_distributed_fft2d_mpi synergia_distributed_fft
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft2d_mpi1 ${MPIEXEC} -np 1 test_distributed_fft2d_mpi)
add_test(test_distributed_fft2d_mpi2 ${MPIEXEC} -np 2 test_distributed_fft2d_mpi)
add_test(test_distributed_fft2d_mpi3 ${MPIEXEC} -np 3 test_distributed_fft2d_mpi)
add_test(test_distributed_fft2d_mpi4 ${MPIEXEC} -np 4 test_distributed_fft2d_mpi)

add_executable(test_xml_serialization test_xml_serialization.cc)
target_link_libraries(test_xml_serialization ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_xml_serialization test_xml_serialization)

add_executable(test_complex_error_function test_complex_error_function.cc)
target_link_libraries(test_complex_error_function ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_complex_error_function test_complex_error_function)

add_executable(test_multi_array_offsets test_multi_array_offsets.cc)
target_link_libraries(test_multi_array_offsets
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_multi_array_offsets test_multi_array_offsets)

add_executable(test_multi_array_serialization test_multi_array_serialization.cc)
target_link_libraries(test_multi_array_serialization ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_multi_array_serialization test_multi_array_serialization)

add_executable(test_command_line_arg test_command_line_arg.cc)
target_link_libraries(test_command_line_arg synergia_command_line
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_command_line_arg test_command_line_arg)

set(CLEAN_FILES
    array1d.h5
    array2d.h5
    array3d.h5
    chunkedarray2d.h5
    double_float.h5
    integer.h5
    justarray1d.h5
    justarray2d.h5
    justarray3d.h5
    justdouble.h5
    justinteger.h5
    hdf5_file_empty.h5
    hdf5_file_write.h5
    a.xml
    multi_array1d.xml
    multi_array2d.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

add_library(wrap_containers wrap_containers.cc)
add_python_extension(wrap_containers_python wrap_containers_wrap.cc)
  target_link_libraries(wrap_containers_python wrap_containers ${Boost_LIBRARIES})
add_test(test_container_conversions.py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_container_conversions.py)

  
#     wrap_containers.h
#      DESTINATION include/synergia/utils/tests)
#   install(TARGETS
#           wrap_containers
#           wrap_containers_python 
#               DESTINATION lib/synergia/utils/tests)


