add_executable(fodo_read_openpmd fodo_cxx.cc fodo_cxx_options.cc)
target_link_libraries(
  fodo_read_openpmd
  synergia_simulation
  synergia_lattice
  synergia_collective
  synergia_serialization
  synergia_bunch
  synergia_command_line)
