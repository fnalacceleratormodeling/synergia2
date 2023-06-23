add_executable(test_libff_elements test_libff_elements.cc)
target_link_libraries(test_libff_elements synergia_simulation
                      synergia_test_main ${kokkos_libs})
add_mpi_test(test_libff_elements 1)
copy_file(fodo.madx test_libff_elements)

add_executable(test_madx_multipoles test_madx_multipoles.cc)
target_link_libraries(test_madx_multipoles synergia_simulation
                      synergia_test_main ${kokkos_libs})
add_mpi_test(test_madx_multipoles 1)

add_executable(test_madx_elements test_madx_elements.cc)
target_link_libraries(test_madx_elements synergia_simulation synergia_test_main
                      ${kokkos_libs})
add_mpi_test(test_madx_elements 1)

add_executable(test_foils test_foils.cc)
target_link_libraries(test_foils synergia_simulation synergia_test_main
                      ${kokkos_libs})
add_mpi_test(test_foils 1)

add_subdirectory(lattices)

add_executable(test_zeroparticle_propagate test_zeroparticle_propagate.cc)
target_link_libraries(test_zeroparticle_propagate synergia_simulation synergia_test_main
                      ${kokkos_libs})
add_mpi_test(test_zeroparticle_propagate 1)
