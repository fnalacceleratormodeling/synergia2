add_custom_target(doxy-bunch doxygen bunch.cfg)
add_dependencies(doxy doxy-bunch)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES xml)
