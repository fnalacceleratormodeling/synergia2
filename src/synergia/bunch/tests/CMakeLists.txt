add_executable(test_bunch_traits test_bunch_traits.cc)
target_link_libraries(test_bunch_traits kokkos)
add_test(NAME test_bunch_traits COMMAND test_bunch_traits)

add_executable(test_bunch test_bunch.cc)
target_link_libraries(test_bunch synergia_bunch)
add_test(NAME test_bunch COMMAND test_bunch)

