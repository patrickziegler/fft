cmake_minimum_required(VERSION 3.12)

project(fft VERSION 0.1.0)

add_subdirectory(src)
add_subdirectory(bindings/matlab)
add_subdirectory(bindings/cython)
add_subdirectory(test/python)
