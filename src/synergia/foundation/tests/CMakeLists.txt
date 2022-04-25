
if (BUILD_PYTHON_BINDINGS)

    add_py_test(test_version.py)
    add_py_test(test_four_momentum.py)
    add_py_test(test_reference_particle.py)

endif()


add_executable(test_trigon test_trigon.cc)
target_link_libraries(test_trigon synergia_foundation ${kokkos_libs} synergia_test_main)
add_mpi_test(test_trigon 1)


add_executable(test_four_momentum test_four_momentum.cc)
target_link_libraries(test_four_momentum synergia_foundation synergia_test_main)
add_mpi_test(test_four_momentum 1)

add_executable(test_reference_particle test_reference_particle.cc)
target_link_libraries(test_reference_particle synergia_foundation synergia_test_main)
add_mpi_test(test_reference_particle 1)
