add_custom_target(update_synergia_cfg
    COMMAND ./update-synergia-cfg)

add_custom_target(doxygen
    /bin/rm -rf xml
    COMMAND doxygen synergia.cfg)
add_dependencies(doxygen update_synergia_cfg)

add_custom_target(update_modules_rst
    COMMAND ./update-modules-rst)

add_custom_target(sphinx
    /bin/rm -rf html doctrees
    COMMAND sphinx-build -a -b html -d doctrees sphinx-src html)
add_dependencies(sphinx doxygen update_modules_rst)

add_custom_target(quickhtml
    /bin/rm -rf html doctrees
    COMMAND sphinx-build -a -b html -d doctrees sphinx-src html)

add_custom_target(sphinx_latex
    /bin/rm -rf pdf doctrees
    COMMAND sphinx-build -a -b latex -d doctrees sphinx-src pdf)
add_dependencies(sphinx_latex doxygen update_modules_rst)

add_custom_target(sphinx_pdf
    COMMAND ./run-pdflatex)
add_dependencies(sphinx_pdf sphinx_latex)

add_custom_target(doc)
add_dependencies(doc sphinx)

set_directory_properties(PROPERTIES ADDITIONAL_MAKE_CLEAN_FILES xml html doctrees)

add_subdirectory(devel/profiling)