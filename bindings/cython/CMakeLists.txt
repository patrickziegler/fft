configure_file(
    setup.py.in
    setup.py
    )

configure_file(
    fft_h.pxd.in
    fft_h.pxd
    )

configure_file(
    fft.pyx
    fft.pyx
    COPYONLY
    )
