add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_four_momentum test_four_momentum.cc)
target_link_libraries(test_four_momentum synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_four_momentum test_four_momentum)
add_test(test_four_momentum_py nosetests test_four_momentum.py)

add_executable(test_reference_particle test_reference_particle.cc)
target_link_libraries(test_reference_particle synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_reference_particle test_reference_particle)
add_test(test_reference_particle_py nosetests test_reference_particle.py)

add_executable(test_distribution test_distribution.cc)
target_link_libraries(test_distribution synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_distribution test_distribution)
add_test(test_distribution_py nosetests test_distribution.py)

set(CLEAN_FILES
    four_momentum.xml
    reference_particle.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

