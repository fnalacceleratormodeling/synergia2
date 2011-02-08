add_definitions(-DBOOST_TEST_DYN_LINK)

add_executable(test_lattice_element test_lattice_element.cc)
target_link_libraries(test_lattice_element synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice_element test_lattice_element)
add_test(test_pylattice_element ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_element.py)

add_executable(test_lattice_element_slice test_lattice_element_slice.cc)
target_link_libraries(test_lattice_element_slice synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice_element_slice test_lattice_element_slice)
# add_test(test_pylattice_element ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_element.py)

add_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice test_lattice)
# add_test(test_pylattice ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice.py)

add_executable(test_chef_lattice test_chef_lattice.cc)
target_link_libraries(test_chef_lattice synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_chef_lattice test_chef_lattice)
# add_test(test_pylattice ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice.py)

add_test(test_mad8_parser ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_parser.py)
add_test(test_expression_parser ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_expression_parser.py)
add_test(test_mad8_reader ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_reader.py)

set(CLEAN_FILES
    lattice_element.xml
    lattice1.xml
    lattice2.xml
    lattice_cache)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")
