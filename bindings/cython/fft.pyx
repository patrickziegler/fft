import numpy as np
cimport numpy as np
cimport cython
cimport fft_h


@cython.boundscheck(False)
@cython.wraparound(False)
def fft_dif_radix2(np.ndarray[np.double_t, ndim=1] x):

    cdef np.ndarray[np.double_t, ndim=1] buf_xr = np.array(x)
    cdef np.ndarray[np.double_t, ndim=1] buf_xi = np.zeros(x.size)

    fft_h.fft_dif_radix2(x.size, &buf_xr[0], &buf_xi[0])

    cdef np.ndarray[np.double_t, ndim=1] buf_Xr = np.zeros(x.size)
    cdef np.ndarray[np.double_t, ndim=1] buf_Xi = np.zeros(x.size)

    fft_h.copy_bitrev(x.size, &buf_xr[0], &buf_xi[0], &buf_Xr[0], &buf_Xi[0])

    return buf_Xr + 1j * buf_Xi


@cython.boundscheck(False)
@cython.wraparound(False)
def fft_dif_radix4(np.ndarray[np.double_t, ndim=1] x):

    cdef np.ndarray[np.double_t, ndim=1] buf_xr = np.array(x)
    cdef np.ndarray[np.double_t, ndim=1] buf_xi = np.zeros(x.size)

    fft_h.fft_dif_radix4(x.size, &buf_xr[0], &buf_xi[0])

    cdef np.ndarray[np.double_t, ndim=1] buf_Xr = np.zeros(x.size)
    cdef np.ndarray[np.double_t, ndim=1] buf_Xi = np.zeros(x.size)

    fft_h.copy_bitrev(x.size, &buf_xr[0], &buf_xi[0], &buf_Xr[0], &buf_Xi[0])

    return buf_Xr + 1j * buf_Xi
