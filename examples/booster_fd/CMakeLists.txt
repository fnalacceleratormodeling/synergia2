add_executable(booster_fd booster_fd.cc)
target_link_libraries(booster_fd
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    synergia_command_line
    )

copy_file(orbump_rampdown_0099.lsx booster_fd)
copy_file(particles_booster_initial_onepercent.h5 booster_fd)
