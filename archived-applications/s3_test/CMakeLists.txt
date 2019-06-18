add_executable(s3_test s3_test.cc)
target_link_libraries(s3_test 
    synergia_simulation 
    synergia_lattice
    synergia_bunch 
    synergia_libFF
    gsl
    MPI::MPI_C
    )

add_custom_command(TARGET s3_test
    POST_BUILD
    COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/fodo.madx 
        ${CMAKE_CURRENT_BINARY_DIR}/fodo.madx
    COMMENT "Copying lattice file for s3_test"
    )


