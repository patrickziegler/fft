#ifndef FFT_H
#define FFT_H

#include <cstddef>

void copy_bitrev(const size_t&, const double*, const double*, double*, double*);
void fft_dif_radix2(const size_t&, double*, double*);
void fft_dif_radix4(const size_t&, double*, double*);

#endif
