add_executable(fodo_cxx fodo_cxx.cc)
target_link_libraries(fodo_cxx
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    )

