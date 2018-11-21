add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_chef_utils test_chef_utils.cc)
target_link_libraries(test_chef_utils synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(test_chef_utils test_chef_utils)

add_test_executable(test_lattice_element test_lattice_element.cc)
target_link_libraries(test_lattice_element synergia_lattice lsexpr
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice_element COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_lattice_element>)
add_test(NAME test_lattice_element_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_element.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_lattice_element.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_lattice_element_slice test_lattice_element_slice.cc)
target_link_libraries(test_lattice_element_slice synergia_lattice
    synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice_element_slice COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_lattice_element_slice>)
# add_test(NAME test_pylattice_element COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_element.py)

add_test_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice synergia_lattice synergia_foundation lsexpr
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_lattice>)
add_test(NAME test_lattice_py COMMAND
    ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1
    ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()"
    test_lattice.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_lattice.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_lattice_diagnostics test_lattice_diagnostics.cc)
target_link_libraries(test_lattice_diagnostics synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice_diagnostics COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_lattice_diagnostics>)

add_test_executable(test_lattice_diagnostics_mpi test_lattice_diagnostics_mpi.cc)
target_link_libraries(test_lattice_diagnostics_mpi synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice_diagnostics_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 $<TARGET_FILE:test_lattice_diagnostics_mpi>)

add_test_executable(test_apertures_diagnostics_mpi3 test_apertures_diagnostics.cc)
target_link_libraries(test_apertures_diagnostics_mpi3 synergia_lattice synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
#add_test(NAME test_apertures_diagnostics COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_apertures_diagnostics>)
add_test(NAME test_apertures_diagnostics_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 $<TARGET_FILE:test_apertures_diagnostics_mpi3>)


add_test_executable(test_element_adaptor test_element_adaptor.cc)
target_link_libraries(test_element_adaptor synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_element_adaptor COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_element_adaptor>)

add_test_executable(test_mad8_adaptor_map test_mad8_adaptor_map.cc)
target_link_libraries(test_mad8_adaptor_map synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_mad8_adaptor_map COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_mad8_adaptor_map>)

add_test_executable(test_madx_adaptor_map test_madx_adaptor_map.cc)
target_link_libraries(test_madx_adaptor_map synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_madx_adaptor_map COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_madx_adaptor_map>)

add_test_executable(test_chef_lattice test_chef_lattice.cc)
target_link_libraries(test_chef_lattice synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_lattice COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_lattice>)
# add_test(NAME test_pylattice COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice.py)

add_test_executable(test_chef_lattice_section test_chef_lattice_section.cc)
target_link_libraries(test_chef_lattice_section synergia_lattice synergia_foundation
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_lattice_section COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_lattice_section>)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/foo.dbx DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test_executable(test_madx_parser test_madx_parser.cc)
target_link_libraries(test_madx_parser synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_madx_parser COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  $<TARGET_FILE:test_madx_parser>)

add_test_executable(test_madx_reader test_madx_reader.cc)
target_link_libraries(test_madx_reader synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_madx_reader COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  $<TARGET_FILE:test_madx_reader>)

add_test_executable(jfa jfa.cc)
target_link_libraries(jfa synergia_lattice
    ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME test_mad8_parser_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_parser.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_mad8_parser.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_expression_parser_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_expression_parser.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_expression_parser.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_mad8_reader_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_reader.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_mad8_reader.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_simplify_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_simplify.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_simplify.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lattices
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
    
set(CLEAN_FILES
    lattice_element.xml
    lattice1.xml
    lattice2.xml
    lattice_cache)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")
