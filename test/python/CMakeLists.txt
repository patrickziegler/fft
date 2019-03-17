find_package(Python3)

configure_file(
    Signal.py
    Signal.py
    COPYONLY
    )

configure_file(
    test_pyfft.py
    test_pyfft.py
    COPYONLY
    )

add_custom_target(
    test_env
    DEPENDS env
    )

add_custom_command(
    OUTPUT env
    COMMAND ${Python3_EXECUTABLE} -m virtualenv env --system-site-packages
    )

add_custom_target(
    test_pyfft ALL
    COMMAND source env/bin/activate && cd "${CMAKE_BINARY_DIR}/bindings/cython" && python setup.py develop --prefix "${CMAKE_CURRENT_BINARY_DIR}/env"
    COMMAND source env/bin/activate && python test_pyfft.py
    )

add_dependencies(test_pyfft test_env)
