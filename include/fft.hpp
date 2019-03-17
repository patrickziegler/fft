#ifndef FFT_H
#define FFT_H

#include <cstddef>

void fft_dif_radix2(const size_t, double*, double*);
void fft_dif_radix4(const size_t, double*, double*);
void copy_bitrev(size_t, double*, double*, double*, double*);

#endif
