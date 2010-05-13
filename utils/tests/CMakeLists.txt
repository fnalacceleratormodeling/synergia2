include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_parallel_utils test_parallel_utils.cc)
target_link_libraries(test_parallel_utils parallel_utils ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_parallel_utils test_parallel_utils)
add_test(test_parallel_utils_py nosetests test_parallel_utils.py)

add_executable(test_floating_point test_floating_point.cc)
target_link_libraries(test_floating_point ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_floating_point test_floating_point)

add_executable(test_hdf5_writer test_hdf5_writer.cc)
target_link_libraries(test_hdf5_writer ${HDF5_LIBRARIES} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_hdf5_writer test_hdf5_writer)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES integer.h5 double_float.h5
    array1d.h5 array2d.h array3d.h)
