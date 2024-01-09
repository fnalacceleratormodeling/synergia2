add_executable(test_bunch_simulator_mpi test_bunch_simulator_mpi.cc)
target_link_libraries(test_bunch_simulator_mpi synergia_simulation
                      synergia_test_main)

add_mpi_test(test_bunch_simulator_mpi 1)
add_mpi_test(test_bunch_simulator_mpi 2)
add_mpi_test(test_bunch_simulator_mpi 3)
add_mpi_test(test_bunch_simulator_mpi 4)
add_mpi_test(test_bunch_simulator_mpi 5)
add_mpi_test(test_bunch_simulator_mpi 6)

add_executable(test_bunch_simulator test_bunch_simulator.cc)
target_link_libraries(test_bunch_simulator synergia_simulation
                      synergia_test_main)
add_mpi_test(test_bunch_simulator 1)

add_executable(test_get_tunes2 test_get_tunes2.cc)
target_link_libraries(test_get_tunes2 synergia_simulation synergia_test_main)
add_mpi_test(test_get_tunes2 1)

add_executable(test_get_tunes3 test_get_tunes3.cc)
target_link_libraries(test_get_tunes3 synergia_simulation synergia_test_main)
add_mpi_test(test_get_tunes3 1)

add_executable(test_get_tunes_dpp_offset test_get_tunes_dpp_offset.cc)
target_link_libraries(test_get_tunes_dpp_offset synergia_simulation
                      synergia_test_main)
add_mpi_test(test_get_tunes_dpp_offset 1)

add_executable(test_closed_orbits_and_others test_closed_orbits_and_others.cc)
target_link_libraries(test_closed_orbits_and_others synergia_simulation
                      synergia_test_main)
add_mpi_test(test_closed_orbits_and_others 1)

add_executable(test_closed_orbit2 test_closed_orbit2.cc)
target_link_libraries(test_closed_orbit2 synergia_simulation
                      synergia_test_main)
add_mpi_test(test_closed_orbit2 1)

add_executable(test_closed_orbit3 test_closed_orbit3.cc)
target_link_libraries(test_closed_orbit3 synergia_simulation
                      synergia_test_main)
add_mpi_test(test_closed_orbit3 1)

if(BUILD_PYTHON_BINDINGS)
  add_py_test(test_propagator.py)
  add_py_test(test_nonlinear_maps.py)
  add_py_test(test_booster_normal_form.py)
  add_py_test(test_rfcavity.py)
  add_py_test(test_rfcavity2.py)
  add_py_test(test_tune_circular_lattice.py)
  add_py_test(test_tune_circular_lattice2.py)
  add_py_test(test_prop_actions.py)
  add_py_test(test_accel.py)
  # test_accel2.py removed as redundant
  add_py_test(test_accel3.py)
  add_py_test(test_accel4.py)
  add_py_test(test_accel5.py)
  add_py_test(test_accel_booster.py)
  add_py_test(test_accel_diagnostics.py)
  add_py_test(test_bend2.py)
  add_py_test(test_higher_momentum_bend.py)
  add_py_test(test_multiple_momentum_propagation.py)
  add_py_test(test_multiple_momentum_propagation_sbend.py)
  add_py_test(test_rfcavity_plus_drift.py)
  add_py_test(test_rfcavity_plus_bend.py)
  add_py_test(test_rfcavity_plus_cfbend.py)
  add_py_test(test_rfcavity_plus_quad.py)
  add_py_test(test_rfcavity_plus_sex.py)
  add_py_test(test_rfcavity_plus_oct.py)
  add_py_test(test_kicked_orbit_timing.py)
  add_py_test(test_kicked_orbit_timing2.py)
  add_py_test(test_get_tunes.py)
  add_py_test(test_empty_matrix_map.py)
  add_py_test(test_matrix_map.py)
  add_py_test(test_matrix_prop.py)
  add_py_test(test_matrix_propkick.py)
endif()
