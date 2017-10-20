add_definitions(-DBOOST_TEST_DYN_LINK)

add_test_executable(test_chef_propagator test_chef_propagator.cc)
target_link_libraries(test_chef_propagator synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_propagator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_propagator>)

add_test_executable(test_fast_mapping_term test_fast_mapping_term.cc)
target_link_libraries(test_fast_mapping_term synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_fast_mapping_term COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_fast_mapping_term>)

add_test_executable(test_fast_mapping test_fast_mapping.cc)
target_link_libraries(test_fast_mapping synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_fast_mapping COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_fast_mapping>)

add_test_executable(test_fast_mapping_operation test_fast_mapping_operation.cc)
target_link_libraries(test_fast_mapping_operation synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_fast_mapping_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_fast_mapping_operation>)

add_test_executable(test_chef_propagate_operation test_chef_propagate_operation.cc)
target_link_libraries(test_chef_propagate_operation synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_propagate_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_propagate_operation>)

add_test_executable(test_chef_map_operation_extractor test_chef_map_operation_extractor.cc)
target_link_libraries(test_chef_map_operation_extractor synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_map_operation_extractor COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_map_operation_extractor>)

add_test_executable(test_chef_propagate_operation_extractor test_chef_propagate_operation_extractor.cc)
target_link_libraries(test_chef_propagate_operation_extractor synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_propagate_operation_extractor COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_propagate_operation_extractor>)

add_test_executable(test_chef_mixed_operation_extractor test_chef_mixed_operation_extractor.cc)
target_link_libraries(test_chef_mixed_operation_extractor synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_chef_mixed_operation_extractor COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_chef_mixed_operation_extractor>)

add_test_executable(test_operation_extractor_map test_operation_extractor_map.cc)
target_link_libraries(test_operation_extractor_map synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_operation_extractor_map COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_operation_extractor_map>)

add_test_executable(test_lattice_simulator test_lattice_simulator.cc)
target_link_libraries(test_lattice_simulator synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice_simulator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_lattice_simulator>)
add_test(NAME test_lattice_simulator_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_lattice_simulator.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_lattice_simulator.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_mad8_multipoles_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_mad8_multipoles.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_mad8_multipoles.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_mad8_misc_elements_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_misc_m8_elements.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_misc_m8_elements.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_madx_misc_elements_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_misc_mx_elements.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_misc_mx_elements.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_madx_multipoles_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_madx_multipoles.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_madx_multipoles.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_quad_mad8_v_madx_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_quad_mad8_v_madx.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_quad_mad8_v_madx.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})


add_test_executable(test_collective_operator test_collective_operator.cc)
target_link_libraries(test_collective_operator synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_collective_operator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_collective_operator>)

add_test_executable(test_independent_operator test_independent_operator.cc)
target_link_libraries(test_independent_operator synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_independent_operator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_independent_operator>)

add_test_executable(test_step test_step.cc)
target_link_libraries(test_step synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_step COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_step>)

add_test_executable(test_split_operator_stepper test_split_operator_stepper.cc)
target_link_libraries(test_split_operator_stepper synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_split_operator_stepper COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_split_operator_stepper>)

add_test_executable(test_independent_stepper test_independent_stepper.cc)
target_link_libraries(test_independent_stepper synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_independent_stepper COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_independent_stepper>)

add_test_executable(test_independent_stepper_elements test_independent_stepper_elements.cc)
target_link_libraries(test_independent_stepper_elements synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_independent_stepper_elements COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_independent_stepper_elements>)

add_test_executable(test_split_operator_stepper_elements test_split_operator_stepper_elements.cc)
target_link_libraries(test_split_operator_stepper_elements synergia_simulation
synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_split_operator_stepper_elements COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_split_operator_stepper_elements>)

add_test_executable(test_split_operator_stepper_choice test_split_operator_stepper_choice.cc)
target_link_libraries(test_split_operator_stepper_choice synergia_simulation
synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_split_operator_stepper_choice COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_split_operator_stepper_choice>)

add_test_executable(test_propagator test_propagator.cc)
target_link_libraries(test_propagator synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_propagator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_propagator>)

add_test_executable(test_resume test_resume.cc)
target_link_libraries(test_resume synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_resume COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_resume>)

add_test_executable(test_lattice_elements_action test_lattice_elements_action.cc)
target_link_libraries(test_lattice_elements_action synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_lattice_elements_action COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_lattice_elements_action>)


add_test_executable(test_finite_aperture_operation test_finite_aperture_operation.cc)
target_link_libraries(test_finite_aperture_operation synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_finite_aperture_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_finite_aperture_operation>)

add_test_executable(test_circular_aperture_operation test_circular_aperture_operation.cc)
target_link_libraries(test_circular_aperture_operation synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_circular_aperture_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_circular_aperture_operation>)

add_test_executable(test_elliptical_aperture_operation test_elliptical_aperture_operation.cc)
target_link_libraries(test_elliptical_aperture_operation synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_elliptical_aperture_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_elliptical_aperture_operation>)
add_test(NAME test_elliptical_aperture_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose ; nose.main()" test_elliptical_apertures.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_elliptical_apertures.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test_executable(test_rectangular_aperture_operation test_rectangular_aperture_operation.cc)
target_link_libraries(test_rectangular_aperture_operation synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_rectangular_aperture_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_rectangular_aperture_operation>)

add_test_executable(test_rectangular_with_ears_aperture_operation test_rectangular_with_ears_aperture_operation.cc)
target_link_libraries(test_rectangular_with_ears_aperture_operation synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_rectangular_with_ears_aperture_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_rectangular_with_ears_aperture_operation>)

add_test_executable(test_polygon_aperture_operation test_polygon_aperture_operation.cc)
target_link_libraries(test_polygon_aperture_operation synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_polygon_aperture_operation COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_polygon_aperture_operation>)

add_test_executable(test_propagator_mpi test_propagator_mpi.cc)
target_link_libraries(test_propagator_mpi synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_propagator_mpi1 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  $<TARGET_FILE:test_propagator_mpi>)
add_test(NAME test_propagator_mpi2 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2  $<TARGET_FILE:test_propagator_mpi>)
add_test(NAME test_propagator_mpi3 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3  $<TARGET_FILE:test_propagator_mpi>)
add_test(NAME test_propagator_mpi4 COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4  $<TARGET_FILE:test_propagator_mpi>)

add_test_executable(test_calculate_closed_orbit test_calculate_closed_orbit.cc)
target_link_libraries(test_calculate_closed_orbit synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_calculate_closed_orbit COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1  $<TARGET_FILE:test_calculate_closed_orbit>)

add_test(NAME test_calculate_closed_orbit_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_calculate_closed_orbit.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_calculate_closed_orbit.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_set_rf_frequency_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_set_rf_frequency.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_set_rf_frequency.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_simplify2_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_simplify2.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_simplify2.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

add_test(NAME test_propagate_resume_mpi_1_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_propagate_resume.py)
file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/test_propagate_resume.py 
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})
add_test(NAME test_propagate_resume_mpi_2_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 2 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_propagate_resume.py)
add_test(NAME test_propagate_resume_mpi_3_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 3 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_propagate_resume.py)
add_test(NAME test_propagate_resume_mpi_4_py COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${SYNERGIA2_BINARY_DIR}/synergia-local -c "import nose; nose.main()" test_propagate_resume.py)

add_test_executable(test_bunch_simulator test_bunch_simulator.cc)
target_link_libraries(test_bunch_simulator synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_bunch_simulator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_bunch_simulator>)

add_test_executable(test_bunch_train_simulator test_bunch_train_simulator.cc)
target_link_libraries(test_bunch_train_simulator synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_bunch_train_simulator COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_bunch_train_simulator>)

add_test_executable(test_beamline_element_dynamics test_beamline_element_dynamics.cc)
target_link_libraries(test_beamline_element_dynamics synergia_lattice synergia_simulation ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_beamline_element_dynamics COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_beamline_element_dynamics>)

add_test_executable(test_acceleration test_acceleration.cc)
target_link_libraries(test_acceleration synergia_simulation synergia_lattice synergia_bunch ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_acceleration COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_acceleration>)

add_test_executable(test_propagate_with_apertures test_propagate_with_apertures.cc)
target_link_libraries(test_propagate_with_apertures synergia_simulation synergia_lattice ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})
add_test(NAME test_propagate_with_apertures COMMAND ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 1 $<TARGET_FILE:test_propagate_with_apertures>)

file(COPY ${CMAKE_CURRENT_SOURCE_DIR}/lattices
    DESTINATION ${CMAKE_CURRENT_BINARY_DIR})

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES
    "test_fast_mapping_term.dat;test_fast_mapping.dat;test_fast_mapping2.dat;test_propagate_per_step.h5;test_propagate_per_turn.h5")
