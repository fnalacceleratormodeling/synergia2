add_executable(s3_sis18 s3_sis18.cc)
target_link_libraries(s3_sis18
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    synergia_libFF
    gsl
    MPI::MPI_C
    )

add_custom_command(OUTPUT sis18.madx
    COMMAND ${CMAKE_COMMAND} -E copy 
        ${CMAKE_CURRENT_SOURCE_DIR}/sis18.madx 
        ${CMAKE_CURRENT_BINARY_DIR}/sis18.madx
    DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/sis18.madx
    COMMENT "Copying lattice file for s3_sis18"
    )

add_custom_target(s3_sis18_lattice ALL
    DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/sis18.madx
    )

add_dependencies(s3_sis18 s3_sis18_lattice)

