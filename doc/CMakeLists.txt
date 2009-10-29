add_custom_target(sphinxdoc-sphinx sphinx-build -a -b html -d .doctrees sphinx-src html)
add_dependencies(sphinxdoc sphinxdoc-sphinx)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES html)
set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES .doctrees)
