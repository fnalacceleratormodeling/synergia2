
if (BUILD_PYTHON_BINDINGS)

    add_py_test(test_lattice.py)

endif()

add_executable(test_lattice test_lattice.cc)
target_link_libraries(test_lattice
    synergia_lattice
    synergia_test_main)
add_mpi_test(test_lattice 1)


add_executable(test_mx_expr test_mx_expr.cc)
target_link_libraries(test_mx_expr
    synergia_lattice
    synergia_test_main
    kokkos)
add_mpi_test(test_mx_expr 1)


add_executable(test_madx_parser test_madx_parser.cc)
target_link_libraries(test_madx_parser
    synergia_lattice
    synergia_test_main
    kokkos)
add_mpi_test(test_madx_parser 1)

copy_file(foo.dbx test_madx_parser)


add_executable(test_dynamic_lattice test_dynamic_lattice.cc)
target_link_libraries(test_dynamic_lattice
    synergia_lattice
    synergia_test_main)
add_mpi_test(test_dynamic_lattice 1)


