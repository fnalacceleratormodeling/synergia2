add_executable(mi mi.cc)
target_link_libraries(mi synergia_simulation synergia_lattice
                      synergia_collective synergia_bunch)

copy_file(mi20_raw.lsx mi)

copy_file(mi.py mi)
copy_file(mi_setup.py mi)
copy_file(mi_fixup.py mi)
copy_file(mi_multibunch_options.py mi)
copy_file(mi20_ra_08182020.lat mi)
