add_executable(s3_test s3_test.cc)
target_link_libraries(s3_test 
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    synergia_libFF
    gsl
    MPI::MPI_C
    MPI::MPI_CXX
    )

add_custom_command(OUTPUT fodo.madx
    COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/fodo.madx 
        ${CMAKE_CURRENT_BINARY_DIR}/fodo.madx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/fodo.madx
    COMMENT "Copying lattice file for s3_test"
    )

add_custom_target(s3_test_lattice ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/fodo.madx
    )

add_dependencies(s3_test s3_test_lattice)

