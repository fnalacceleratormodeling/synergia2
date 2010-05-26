# - FindCHEF
# Find CHEF includes and library
# This module defines:
# CHEF_INCLUDE_DIRS
# CHEF_LIBARARY_DIRS
# CHEF_LIBS
# CHEF_FOUND

if(NOT CHEF_FOUND)
    set(INTERNAL_CHEF_FOUND TRUE)
    execute_process(COMMAND "chef-config.sh" "--includes_list"
        OUTPUT_VARIABLE CHEF_INCLUDE_DIRS
        RESULT_VARIABLE CHEF_CONFIG_FAILED
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CHEF_CONFIG_FAILED)
        set(INTERNAL_CHEF_FOUND FALSE)
    else(CHEF_CONFIG_FAILED)
        set (CHEF_INCLUDE_DIRS ${CHEF_INCLUDE_DIRS} CACHE FILEPATH "CHEF include paths")
    endif(CHEF_CONFIG_FAILED)

    execute_process(COMMAND "chef-config.sh" "--lib_dirs_list"
        OUTPUT_VARIABLE CHEF_LIBRARY_DIRS
        RESULT_VARIABLE CHEF_CONFIG_FAILED
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CHEF_CONFIG_FAILED)
        set(INTERNAL_CHEF_FOUND FALSE)
    else(CHEF_CONFIG_FAILED)
        set(CHEF_LIBRARY_DIRS ${CHEF_LIBRARY_DIRS} CACHE FILEPATH "CHEF library paths")
    endif(CHEF_CONFIG_FAILED)

    execute_process(COMMAND "chef-config.sh" "--libs_list"
        OUTPUT_VARIABLE CHEF_LIBS
        RESULT_VARIABLE CHEF_CONFIG_FAILED
        OUTPUT_STRIP_TRAILING_WHITESPACE)
    if(CHEF_CONFIG_FAILED)
        set(INTERNAL_CHEF_FOUND FALSE)
    else(CHEF_CONFIG_FAILED)
        set(CHEF_LIBS ${CHEF_LIBS} CACHE FILEPATH "CHEF libraries")
    endif(CHEF_CONFIG_FAILED)
    set(CHEF_FOUND ${INTERNAL_CHEF_FOUND} CACHE BOOL "CHEF configuration successful" FORCE)
endif(NOT CHEF_FOUND)

include(FindPackageHandleStandardArgs)
find_package_handle_standard_args(CHEF DEFAULT_MSG CHEF_FOUND CHEF_INCLUDE_DIRS CHEF_LIBRARY_DIRS CHEF_LIBS)

