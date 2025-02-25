add_custom_target(cold_fodo_simple)

copy_file(cold_fodo.py cold_fodo_simple)
copy_file (cold_fodo_options.py cold_fodo_simple)
copy_file (fodo.lat cold_fodo_simple)

set(CLEAN_FILES
    log)

set_directory_properties(PROPERTIES ADDITIONAL_CLEAN_FILES "${CLEAN_FILES}")
