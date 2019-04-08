find_package(OpenMP REQUIRED)

add_library(fft SHARED)

target_include_directories(fft PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_sources(fft PRIVATE
    fft.cpp
    )

target_link_libraries(fft OpenMP::OpenMP_CXX)

set_property(TARGET fft PROPERTY POSITION_INDEPENDENT_CODE ON)
