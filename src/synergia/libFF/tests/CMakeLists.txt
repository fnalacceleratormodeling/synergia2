
add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_beamline_element_dynamics_libff test_beamline_element_dynamics_libff.cc)
target_link_libraries(test_beamline_element_dynamics_libff synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_beamline_element_dynamics_libff COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_beamline_element_dynamics_libff>)


add_test_executable(test_libff_drift test_libff_drift.cc)
target_link_libraries(test_libff_drift synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_drift COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_drift>)


