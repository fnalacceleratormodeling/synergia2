if(NOT PYTHON_INCLUDE_DIR)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
        "-c" "import distutils.sysconfig; print distutils.sysconfig.get_python_inc()"
        OUTPUT_VARIABLE PYTHON_INCLUDE_DIR
        RESULT_VARIABLE PYTHON_INCLUDE_NOT_FOUND
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(PYTHON_INCLUDE_NOT_FOUND)
        set(INTERNAL_PYTHON_INCLUDE_FOUND FALSE)
    else(PYTHON_INCLUDE_NOT_FOUND)
        if(PYTHON_INCLUDE_DIR MATCHES "Traceback")
            # Did not successfully find PYTHON_INCLUDE
           set(INTERNAL_PYTHON_INCLUDE_FOUND FALSE)
        else(PYTHON_INCLUDE_DIR MATCHES "Traceback")
            # successful
          set(INTERNAL_PYTHON_INCLUDE_FOUND TRUE)
          set(PYTHON_INCLUDE_DIR ${PYTHON_INCLUDE_DIR} CACHE STRING "Python include path")
        endif (PYTHON_INCLUDE_DIR MATCHES "Traceback")
    endif(PYTHON_INCLUDE_NOT_FOUND)
    set(PYTHON_INCLUDE_FOUND ${INTERNAL_PYTHON_INCLUDE_FOUND} CACHE BOOL "PYTHON_INCLUDE found" FORCE)
endif(NOT PYTHON_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PythonInclude DEFAULT_MSG PYTHON_INCLUDE_DIR)
