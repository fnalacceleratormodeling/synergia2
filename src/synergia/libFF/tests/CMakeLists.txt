
add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_beamline_element_dynamics_libff test_beamline_element_dynamics_libff.cc)
target_link_libraries(test_beamline_element_dynamics_libff synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_beamline_element_dynamics_libff COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_beamline_element_dynamics_libff>)


add_test_executable(test_libff_drift test_libff_drift.cc)
target_link_libraries(test_libff_drift synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_drift COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_drift>)

add_test_executable(test_libff_mem_align test_libff_mem_align.cc)
target_link_libraries(test_libff_mem_align synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_mem_align COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_mem_align>)

add_test_executable(test_libff_on_chef test_libff_on_chef.cc)
target_link_libraries(test_libff_on_chef synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_on_chef COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_on_chef>)

add_test_executable(test_libff_on_iota test_libff_on_iota.cc)
target_link_libraries(test_libff_on_iota synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_on_iota COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_on_iota>)

add_test_executable(test_libff_kicker_on_chef test_libff_kicker_on_chef.cc)
target_link_libraries(test_libff_kicker_on_chef synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_kicker_on_chef COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_kicker_on_chef>)

add_test_executable(test_libff_rfcavities test_libff_rfcavities.cc)
target_link_libraries(test_libff_rfcavities synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_rfcavities COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_rfcavities>)

add_test_executable(test_sis18 test_sis18.cc)
target_link_libraries(test_sis18 synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_sis18 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_sis18>)

add_test_executable(test_libff_acceleration test_libff_acceleration.cc)
target_link_libraries(test_libff_acceleration synergia_lattice synergia_simulation synergia_libFF ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_libff_acceleration COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_libff_acceleration>)

add_test(NAME test_libff_mad8_multipoles_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_multipoles.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_mad8_multipoles.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_libff_madx_multipoles_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_madx_multipoles.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_madx_multipoles.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_nllensff_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_nllens.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_nllens.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
