add_library(fft SHARED)

target_include_directories(fft PUBLIC ${CMAKE_SOURCE_DIR}/include)

target_sources(fft PRIVATE
    fft.cpp
    )

set_property(TARGET fft PROPERTY POSITION_INDEPENDENT_CODE ON)
