add_executable(normal_form normal_form.cc)
target_link_libraries(normal_form
    synergia_simulation 
    synergia_lattice
    synergia_bunch 
    synergia_pydiag
    )

copy_file(channel.madx normal_form)
copy_file(normal_form.py normal_form)

