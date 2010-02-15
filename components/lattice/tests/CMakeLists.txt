add_definitions(-DBOOST_TEST_DYN_LINK)

link_directories(..)

add_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice -llattice -lboost_unit_test_framework)
add_test(test_lattice test_lattice)
add_test(test_pylattice nosetests test_lattice.py)

add_test(test_mad8_parser nosetests test_mad8_parser.py)
add_test(test_expression_parser nosetests test_expression_parser.py)
