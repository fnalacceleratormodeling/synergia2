add_definitions(-DBOOST_TEST_DYN_LINK)

link_directories(..)
link_directories(../../foundation)

add_executable(test_lattice_element test_lattice_element.cc)
target_link_libraries(test_lattice_element -llattice -lboost_unit_test_framework)
add_test(test_lattice_element test_lattice_element)
add_test(test_pylattice_element nosetests test_lattice_element.py)

add_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice -llattice -lfoundation -lboost_unit_test_framework)
# add_test(test_lattice test_lattice)
# add_test(test_pylattice nosetests test_lattice.py)

add_test(test_mad8_parser nosetests test_mad8_parser.py)
add_test(test_expression_parser nosetests test_expression_parser.py)
