include_directories(${Boost_INCLUDES})
add_definitions(-DBOOST_TEST_DYN_LINK)

link_directories(.. ../../foundation)

add_executable(test_bunch test_bunch.cc)
target_link_libraries(test_bunch -lbunch -lfoundation -lboost_unit_test_framework)
add_test(test_bunch test_bunch)
add_test(test_bunch_py nosetests test_bunch.py)

add_executable(test_diagnostics test_diagnostics.cc)
target_link_libraries(test_diagnostics -lbunch -lfoundation -lboost_unit_test_framework)
add_test(test_diagnostics test_diagnostics)
