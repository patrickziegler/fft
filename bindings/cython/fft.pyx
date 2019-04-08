# Copyright (C) 2019 Patrick Ziegler
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.


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
