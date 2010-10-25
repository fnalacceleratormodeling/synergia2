include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_parallel_utils test_parallel_utils.cc)
target_link_libraries(test_parallel_utils synergia_parallel_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_parallel_utils test_parallel_utils)
add_test(test_parallel_utils_py nosetests test_parallel_utils.py)

add_executable(test_floating_point test_floating_point.cc)
target_link_libraries(test_floating_point ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_floating_point test_floating_point)

add_executable(test_hdf5_writer test_hdf5_writer.cc)
target_link_libraries(test_hdf5_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_writer test_hdf5_writer)

add_executable(test_hdf5_serial_writer test_hdf5_serial_writer.cc)
target_link_libraries(test_hdf5_serial_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_serial_writer test_hdf5_serial_writer)

add_executable(test_hdf5_chunked_array2d_writer test_hdf5_chunked_array2d_writer.cc)
target_link_libraries(test_hdf5_chunked_array2d_writer synergia_hdf5_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_chunked_array2d_writer test_hdf5_chunked_array2d_writer)

add_executable(test_distributed_fft3d test_distributed_fft3d.cc)
target_link_libraries(test_distributed_fft3d synergia_distributed_fft ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft3d test_distributed_fft3d)

add_executable(test_distributed_fft3d_mpi test_distributed_fft3d_mpi.cc)
target_link_libraries(test_distributed_fft3d_mpi synergia_distributed_fft ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distributed_fft3d_mpi ${MPIEXEC} -np 2 test_distributed_fft3d_mpi)

add_executable(test_xml_serialization test_xml_serialization.cc)
target_link_libraries(test_xml_serialization ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_xml_serialization test_xml_serialization)

add_executable(test_multi_array_serialization test_multi_array_serialization.cc)
target_link_libraries(test_multi_array_serialization ${Boost_SERIALIZATION_LIBRARY}
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_multi_array_serialization test_multi_array_serialization)

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
    a.xml
    multi_array1d.xml
    multi_array2d.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

