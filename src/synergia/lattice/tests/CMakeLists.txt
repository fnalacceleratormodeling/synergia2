add_executable(test_madx_parser test_madx_parser.cc)
target_link_libraries(test_madx_parser
    synergia_lattice
    synergia_test_main
    kokkos)
add_mpi_test(test_madx_parser 1)

copy_file(foo.dbx test_madx_parser)


