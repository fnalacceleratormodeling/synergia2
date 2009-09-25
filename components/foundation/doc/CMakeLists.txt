add_custom_target(doxy-foundation doxygen foundation.cfg)
add_dependencies(doxy doxy-foundation)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES xml)
