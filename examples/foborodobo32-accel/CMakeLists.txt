add_executable(pltest pltest.cc)
target_link_libraries(pltest
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    synergia_command_line
    )

