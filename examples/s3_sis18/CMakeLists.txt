add_executable(s3_sis18 s3_sis18.cc)
target_link_libraries(s3_sis18
    synergia_simulation 
    synergia_lattice
    synergia_collective
    synergia_bunch 
    synergia_libFF
    )

copy_file(sis18.madx s3_sis18)
copy_file(sis18-6.lsx s3_sis18)
copy_file(sis18.py s3_sis18)
copy_file(resume.py s3_sis18)

