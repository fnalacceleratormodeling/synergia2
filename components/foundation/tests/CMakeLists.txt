include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

include_directories(..)
link_directories(..)

add_executable(test_four_momentum test_four_momentum.cc)
target_link_libraries(test_four_momentum -lfoundation -lboost_unit_test_framework)
add_test(test_four_momentum test_four_momentum)
add_test(test_four_momentum_py nosetests test_four_momentum.py)

add_executable(test_reference_particle test_reference_particle.cc)
target_link_libraries(test_reference_particle -lfoundation -lboost_unit_test_framework)
add_test(test_reference_particle test_reference_particle)
add_test(test_reference_particle_py nosetests test_reference_particle.py)

add_executable(test_distribution test_distribution.cc)
target_link_libraries(test_distribution -lfoundation -lboost_unit_test_framework)
add_test(test_distribution test_distribution)
# add_test(test_distribution_py nosetests test_distribution.py)
