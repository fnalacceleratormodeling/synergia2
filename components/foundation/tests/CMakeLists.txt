include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

#add_executable(unit_test_example_01 unit_test_example_01.cpp)
#target_link_libraries(unit_test_example_01 -lboost_unit_test_framework)
#add_test(unit_test_example_01 unit_test_example_01)

include_directories(..)
link_directories(..)

add_executable(test_four_momentum test_four_momentum.cc)
target_link_libraries(test_four_momentum -lfoundation -lboost_unit_test_framework)
add_test(test_four_momentum test_four_momentum)
 