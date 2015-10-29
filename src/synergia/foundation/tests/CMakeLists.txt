add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_four_momentum test_four_momentum.cc)
target_link_libraries(test_four_momentum synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_four_momentum COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_four_momentum>)
add_test(NAME test_four_momentum_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_four_momentum.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_four_momentum.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_physical_constants test_physical_constants.cc)
target_link_libraries(test_physical_constants synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_physical_constants COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_physical_constants>)

add_test_executable(test_reference_particle test_reference_particle.cc)
target_link_libraries(test_reference_particle synergia_foundation lsexpr ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_reference_particle COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_reference_particle>)
add_test(NAME test_reference_particle_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_reference_particle.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_reference_particle.py
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_distribution test_distribution.cc)
target_link_libraries(test_distribution synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_distribution COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_distribution>)
add_test(NAME test_distribution_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_distribution.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_distribution.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_diagnostics_write_helper test_diagnostics_write_helper.cc)
target_link_libraries(test_diagnostics_write_helper synergia_hdf5_utils synergia_foundation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_diagnostics_write_helper COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_diagnostics_write_helper>)

set(CLEAN_FILES
    four_momentum.xml
    reference_particle.xml)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES "${CLEAN_FILES}")

