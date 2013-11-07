add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_commxx_mpi test_commxx_mpi.cc)
target_link_libraries(test_commxx_mpi synergia_parallel_utils synergia_serialization ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_commxx_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_commxx_mpi)
add_test(test_commxx_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_commxx_mpi)
add_test(test_commxx_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_commxx_mpi)
add_test(test_commxx_mpi7 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 7 test_commxx_mpi)
add_test(test_commxx_mpi12 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 12 test_commxx_mpi)

add_test_executable(test_commxx_divider_mpi test_commxx_divider_mpi.cc)
target_link_libraries(test_commxx_divider_mpi synergia_parallel_utils synergia_serialization ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_commxx_divider_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_commxx_divider_mpi)
add_test(test_commxx_divider_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_commxx_divider_mpi)


add_test_executable(test_parallel_utils test_parallel_utils.cc)
target_link_libraries(test_parallel_utils synergia_parallel_utils
     ${Boost_SERIALIZATION_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_parallel_utils ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_parallel_utils)
add_test(test_parallel_utils_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()"
    test_parallel_utils.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_parallel_utils.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_logger test_logger.cc)
target_link_libraries(test_logger synergia_parallel_utils  ${Boost_SERIALIZATION_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_logger ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_logger)
add_test(test_logger_py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_logger.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_logger.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_logger_mpi test_logger_mpi.cc)
target_link_libraries(test_logger_mpi synergia_parallel_utils ${Boost_SERIALIZATION_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_logger_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_logger_mpi)
add_test(test_logger_mpi3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 test_logger_mpi)

add_test_executable(test_floating_point test_floating_point.cc)
target_link_libraries(test_floating_point ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_floating_point test_floating_point)

add_test_executable(test_fast_int_floor test_fast_int_floor.cc)
target_link_libraries(test_fast_int_floor ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_fast_int_floor test_fast_int_floor)

add_test_executable(test_hdf5_writer test_hdf5_writer.cc)
target_link_libraries(test_hdf5_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_writer test_hdf5_writer)

add_test_executable(test_hdf5_file test_hdf5_file.cc)
target_link_libraries(test_hdf5_file synergia_hdf5_utils synergia_serialization ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_file ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_hdf5_file)

add_test_executable(test_hdf5_serial_writer test_hdf5_serial_writer.cc)
target_link_libraries(test_hdf5_serial_writer synergia_hdf5_utils synergia_serialization ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_serial_writer ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_hdf5_serial_writer)

add_test_executable(test_hdf5_chunked_array2d_writer test_hdf5_chunked_array2d_writer.cc)
target_link_libraries(test_hdf5_chunked_array2d_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_chunked_array2d_writer test_hdf5_chunked_array2d_writer)

add_test_executable(test_distributed_fft3d test_distributed_fft3d.cc)
target_link_libraries(test_distributed_fft3d synergia_distributed_fft ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft3d ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_distributed_fft3d)

add_test_executable(test_distributed_fft3d_mpi test_distributed_fft3d_mpi.cc)
target_link_libraries(test_distributed_fft3d_mpi synergia_distributed_fft
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft3d_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_distributed_fft3d_mpi)
add_test(test_distributed_fft3d_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_distributed_fft3d_mpi)
add_test(test_distributed_fft3d_mpi3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 test_distributed_fft3d_mpi)
add_test(test_distributed_fft3d_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_distributed_fft3d_mpi)

add_test_executable(test_distributed_fft2d test_distributed_fft2d.cc)
target_link_libraries(test_distributed_fft2d synergia_distributed_fft ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft2d ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_distributed_fft2d)

add_test_executable(test_distributed_fft2d_mpi test_distributed_fft2d_mpi.cc)
target_link_libraries(test_distributed_fft2d_mpi synergia_distributed_fft
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft2d_mpi1 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_distributed_fft2d_mpi)
add_test(test_distributed_fft2d_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_distributed_fft2d_mpi)
add_test(test_distributed_fft2d_mpi3 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 test_distributed_fft2d_mpi)
add_test(test_distributed_fft2d_mpi4 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 test_distributed_fft2d_mpi)

add_test_executable(test_xml_serialization test_xml_serialization.cc)
target_link_libraries(test_xml_serialization synergia_serialization
    ${Boost_SERIALIZATION_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_xml_serialization ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_xml_serialization)

add_test_executable(test_complex_error_function test_complex_error_function.cc)
target_link_libraries(test_complex_error_function ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_complex_error_function test_complex_error_function)

add_test_executable(test_multi_array_offsets test_multi_array_offsets.cc)
target_link_libraries(test_multi_array_offsets
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_multi_array_offsets test_multi_array_offsets)

add_test_executable(test_multi_array_serialization test_multi_array_serialization.cc)
target_link_libraries(test_multi_array_serialization synergia_serialization
    ${Boost_SERIALIZATION_LIBRARY} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_multi_array_serialization ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_multi_array_serialization)

add_test_executable(test_command_line_arg test_command_line_arg.cc)
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

# add_library(wrap_containers wrap_containers.cc)
# add_python_extension(wrap_containers_python wrap_containers_wrap.cc)
#   target_link_libraries(wrap_containers_python wrap_containers ${Boost_LIBRARIES})
# add_test(test_container_conversions.py ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_container_conversions.py)


#     wrap_containers.h
#      DESTINATION include/synergia/utils/tests)
#   install(TARGETS
#           wrap_containers
#           wrap_containers_python
#               DESTINATION lib/synergia/utils/tests)


