add_executable(allreduce_benchmark allreduce_benchmark.cc
                                   allreduce_benchmark_options.cc)
target_link_libraries(allreduce_benchmark synergia_parallel_utils
                      synergia_serialization synergia_command_line)
