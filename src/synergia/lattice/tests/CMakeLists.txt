add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_chef_utils test_chef_utils.cc)
target_link_libraries(test_chef_utils synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_chef_utils test_chef_utils)

add_test_executable(test_lattice_element test_lattice_element.cc)
target_link_libraries(test_lattice_element synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice_element ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_lattice_element)
add_test(test_pylattice_element ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_element.py)

add_test_executable(test_lattice_element_slice test_lattice_element_slice.cc)
target_link_libraries(test_lattice_element_slice synergia_lattice
    synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice_element_slice ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_lattice_element_slice)
# add_test(test_pylattice_element ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_element.py)

add_test_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_lattice)
# add_test(test_pylattice ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice.py)

add_test_executable(test_lattice_diagnostics test_lattice_diagnostics.cc)
target_link_libraries(test_lattice_diagnostics synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice_diagnostics ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_lattice_diagnostics)

add_test_executable(test_lattice_diagnostics_mpi test_lattice_diagnostics_mpi.cc)
target_link_libraries(test_lattice_diagnostics_mpi synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_lattice_diagnostics_mpi2 ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 test_lattice_diagnostics_mpi)

add_test_executable(test_element_adaptor test_element_adaptor.cc)
target_link_libraries(test_element_adaptor synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_element_adaptor ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_element_adaptor)

add_test_executable(test_mad8_adaptor_map test_mad8_adaptor_map.cc)
target_link_libraries(test_mad8_adaptor_map synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_mad8_adaptor_map ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_mad8_adaptor_map)

add_test_executable(test_madx_adaptor_map test_madx_adaptor_map.cc)
target_link_libraries(test_madx_adaptor_map synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_madx_adaptor_map ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_madx_adaptor_map)

add_test_executable(test_chef_lattice test_chef_lattice.cc)
target_link_libraries(test_chef_lattice synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_chef_lattice ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_chef_lattice)
# add_test(test_pylattice ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice.py)

add_test_executable(test_chef_lattice_section test_chef_lattice_section.cc)
target_link_libraries(test_chef_lattice_section synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_chef_lattice_section ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 test_chef_lattice_section)

add_test_executable(test_madx_parser test_madx_parser.cc)
target_link_libraries(test_madx_parser synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_madx_parser ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  test_madx_parser)

add_test_executable(test_madx_reader test_madx_reader.cc)
target_link_libraries(test_madx_reader synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_madx_reader ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  test_madx_reader)

add_test_executable(jfa jfa.cc)
target_link_libraries(jfa synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(test_mad8_parser ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_parser.py)
add_test(test_expression_parser ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_expression_parser.py)
add_test(test_mad8_reader ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_reader.py)
add_test(test_simplify ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_simplify.py)

set(CLEAN_FILES
    lattice_element.xml
    lattice1.xml
    lattice2.xml
    lattice_cache)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")
