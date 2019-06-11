add_executable(s3_test s3_test.cc)
target_link_libraries(s3_test 
    synergia_simulation 
    synergia_lattice
    synergia_bunch 
    synergia_libFF
    gsl
    MPI::MPI_C
    )

