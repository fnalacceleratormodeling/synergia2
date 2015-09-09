if(NOT PYTHON_LIB_DIR)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
        "-c" "import distutils.sysconfig; import os.path; print os.path.dirname(distutils.sysconfig.get_python_lib(standard_lib=1))"
        OUTPUT_VARIABLE PYTHON_LIB_DIR
        RESULT_VARIABLE PYTHON_LIB_NOT_FOUND
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(PYTHON_LIB_NOT_FOUND)
        set(INTERNAL_PYTHON_LIB_FOUND FALSE)
    else(PYTHON_LIB_NOT_FOUND)
        if(PYTHON_LIB_DIR MATCHES "Traceback")
            # Did not successfully find PYTHON_LIB
           set(INTERNAL_PYTHON_LIB_FOUND FALSE)
        else(PYTHON_LIB_DIR MATCHES "Traceback")
            # successful
          set(INTERNAL_PYTHON_LIB_FOUND TRUE)
          set(PYTHON_LIB_DIR ${PYTHON_LIB_DIR} CACHE STRING "Python LIB path")
        endif (PYTHON_LIB_DIR MATCHES "Traceback")
    endif(PYTHON_LIB_NOT_FOUND)
    set(PYTHON_LIB_FOUND ${INTERNAL_PYTHON_LIB_FOUND} CACHE BOOL "PYTHON_LIB found" FORCE)
endif(NOT PYTHON_LIB_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(PythonLibDir DEFAULT_MSG PYTHON_LIB_DIR)
