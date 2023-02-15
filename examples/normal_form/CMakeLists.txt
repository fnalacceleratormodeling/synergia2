add_executable(normal_form normal_form.cc)
target_link_libraries(normal_form synergia_simulation synergia_lattice
                      synergia_bunch)

copy_file(channel.madx normal_form)
copy_file(normal_form.py normal_form)

add_executable(channel_map channel_map.cc)
target_link_libraries(channel_map synergia_lattice synergia_simulation)

copy_file(channel_map.py channel_map)
copy_file(foborodobo128.madx channel_map)
