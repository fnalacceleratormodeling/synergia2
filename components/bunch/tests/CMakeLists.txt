include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_bunch test_bunch.cc)
target_link_libraries(test_bunch bunch foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_bunch test_bunch)
add_test(test_bunch_py nosetests test_bunch.py)

add_custom_command(OUTPUT test_diagnostics_get_mean.icc
    COMMAND python test_diagnostics_crosscheck.py
    DEPENDS test_diagnostics_crosscheck.py)
add_executable(test_diagnostics test_diagnostics.cc)
add_dependencies(test_diagnostics test_diagnostics_get_mean.icc)
target_link_libraries(test_diagnostics bunch foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_diagnostics test_diagnostics)
add_test(test_diagnostics_py nosetests test_diagnostics.py)

add_executable(test_populate test_populate.cc)
target_link_libraries(test_populate bunch foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_populate test_populate)
#add_test(test_populate_py nosetests test_populate.py)
