add_custom_target(doxy-utils doxygen utils.cfg)
add_dependencies(doxy doxy-utils)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES xml)
