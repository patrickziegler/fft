set(Matlab_ROOT_DIR "~/.opt/MATLAB/R2017a")

find_package(Matlab)

if(Matlab_FOUND)

    matlab_add_mex(
        NAME fft_dif_radix2_mex
        SRC fft_dif_radix2_mex.cpp
        LINK_TO fft
        MODULE
        )

    matlab_add_mex(
        NAME fft_dif_radix4_mex
        SRC fft_dif_radix4_mex.cpp
        LINK_TO fft
        MODULE
        )

endif()
