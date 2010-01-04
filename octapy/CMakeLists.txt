
add_custom_command(OUTPUT octave_class.o
    DEPENDS octave_class.cc octave_class.h
    COMMAND mkoctfile -I${PYTHON_INCLUDE_PATH} -c octave_class.cc)

add_custom_command(OUTPUT octapy_wrap.cc
    DEPENDS octapy.i octave_class.h
    COMMAND swig -python -c++ -o octapy_wrap.cc octapy.i)

add_custom_command(OUTPUT octapy_wrap.o
    DEPENDS octapy_wrap.cc
    COMMAND mkoctfile -I${PYTHON_INCLUDE_PATH} -c octapy_wrap.cc)

add_custom_command(OUTPUT _octapy.so
    DEPENDS octave_class.o octapy_wrap.o
    COMMAND mkoctfile -I${PYTHON_INCLUDE_PATH} octave_class.o octapy_wrap.o -o _octapy.so)

add_custom_target(_octapy DEPENDS _octapy.so)
