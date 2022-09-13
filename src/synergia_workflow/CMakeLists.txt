add_custom_target(workflow)

copy_file(__init__.py workflow)
copy_file(options.py workflow)
copy_file(job_manager.py workflow)

install(FILES __init__.py options.py job_manager.py
        DESTINATION ${PYTHON_INSTALL_DIR}/synergia_workflow)
