add_executable(benchmark benchmark.cc benchmark_options.cc)
target_link_libraries(benchmark synergia_simulation synergia_collective synergia_command_line)
copy_file(fodo.madx benchmark)
copy_file(cxx_means.txt benchmark)
copy_file(cxx_covariance_matrix.txt benchmark)

add_executable(benchmark2d benchmark2d.cc benchmark_options.cc)
target_link_libraries(benchmark2d synergia_simulation synergia_collective synergia_command_line)

add_executable(sis18 s3_sis18.cc)
target_link_libraries(sis18 synergia_simulation synergia_collective)

copy_file(sis18-6.lsx sis18)

add_executable(read_xml read_xml.cc)
target_link_libraries(read_xml synergia_simulation)
copy_file(read_xml benchmark)
copy_file(write_xml.py benchmark)

add_executable(cxx_write_xml cxx_write_xml.cc)
target_link_libraries(cxx_write_xml synergia_simulation)
copy_file(cxx_write_xml benchmark)
