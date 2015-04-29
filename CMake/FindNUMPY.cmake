# - FindNUMPY
# Find Numpy includes and library
# This module defines:
# NUMPY_INCLUDE_DIR
# NUMPY_FOUND

if(NOT NUMPY_INCLUDE_DIR)
    if(EXTRA_PYTHONPATH)
        set(OLD_PYTHONPATH $ENV{PYTHONPATH})
        set(ENV{PYTHONPATH} "${EXTRA_PYTHONPATH}:$ENV{PYTHONPATH}")
    endif(EXTRA_PYTHONPATH)
    execute_process(COMMAND "${PYTHON_EXECUTABLE}"
        "-c" "import numpy; print numpy.get_include()"
        OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
        RESULT_VARIABLE NUMPY_NOT_FOUND
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(NUMPY_NOT_FOUND)
        set(INTERNAL_NUMPY_FOUND FALSE)
    else(NUMPY_NOT_FOUND)
        if(NUMPY_INCLUDE_DIR MATCHES "Traceback")
            # Did not successfully include numpy
           set(INTERNAL_NUMPY_FOUND FALSE)
        else(NUMPY_INCLUDE_DIR MATCHES "Traceback")
            # successful
          set(INTERNAL_NUMPY_FOUND TRUE)
          set(NUMPY_INCLUDE_DIR ${NUMPY_INCLUDE_DIR} CACHE STRING "Numpy include path")
        endif (NUMPY_INCLUDE_DIR MATCHES "Traceback")
    endif(NUMPY_NOT_FOUND)
    set(NUMPY_FOUND ${INTERNAL_NUMPY_FOUND} CACHE BOOL "Numpy found" FORCE)
    if(EXTRA_PYTHONPATH)
        set(ENV{PYTHONPATH} ${OLD_PYTHONPATH})
    endif(EXTRA_PYTHONPATH)
endif(NOT NUMPY_INCLUDE_DIR)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(NUMPY DEFAULT_MSG NUMPY_INCLUDE_DIR)


