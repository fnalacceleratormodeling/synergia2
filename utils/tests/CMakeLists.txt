include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

include_directories(..)
link_directories(..)

add_executable(test_parallel_utils test_parallel_utils.cc)
target_link_libraries(test_parallel_utils -lparallel_utils -lboost_unit_test_framework)
add_test(test_parallel_utils test_parallel_utils)
add_test(test_parallel_utils_py nosetests test_parallel_utils.py)

add_executable(test_floating_point test_floating_point.cc)
target_link_libraries(test_floating_point -lboost_unit_test_framework)
add_test(test_floating_point test_floating_point)
