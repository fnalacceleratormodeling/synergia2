macro(copy_file fname target)

    #message(${fname})
    #message(${target})

    add_custom_command(OUTPUT ${fname}
        COMMAND ${CMAKE_COMMAND} -E copy
            ${CMAKE_CURRENT_SOURCE_DIR}/${fname}
            ${CMAKE_CURRENT_BINARY_DIR}/${fname}
        DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/${fname}
        COMMENT "Copying file ${fname}"
        )

    add_custom_target(${target}_${fname} ALL
        DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/${fname}
        )

    add_dependencies(${target} ${target}_${fname})

endmacro()



