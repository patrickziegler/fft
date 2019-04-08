// Copyright (C) 2019 Patrick Ziegler
//
// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <https://www.gnu.org/licenses/>.


#include "fft.hpp"
#include <complex>
#include <tuple>

void copy_bitrev(const size_t& n, const double* in_r, const double* in_i, double* out_r, double* out_i)
{
    size_t tmp;
    size_t m;
    size_t n2 = n / 2;
    size_t* idx = new size_t[n];
    idx[0] = 0;

    for (size_t i = 1; i < n; ++i)
    {
        m = n2;
        tmp = idx[i-1];

        while (tmp & m)
        {
            tmp ^= m;
            m >>= 1;
        }

        idx[i] = tmp | m;
    }

    for (size_t i = 0; i < n; ++i)
    {
        out_r[i] = in_r[idx[i]];
        out_i[i] = in_i[idx[i]];
    }

    delete[] idx;
}

std::tuple<size_t, size_t*, double**, double**> complex_factors(const size_t& n, const size_t& radix)
{
    size_t m1;
    size_t m2 = n;

    std::complex<double> e0;
    std::complex<double> e1;

    size_t M = static_cast<size_t>(log2(n) / log2(radix));
    size_t *m = new size_t[M];
    double **er = new double*[M];
    double **ei = new double*[M];

    for (size_t i = 0; i < M; ++i)
    {
        m1 = m2;
        m2 /= radix;
        m[i] = m2;

        er[i] = new double[m2];
        ei[i] = new double[m2];

        e0 = std::exp(std::complex<double>(0, -2 * M_PI / m1));
        e1 = 1;

        for (size_t j = 0; j < m2; ++j)
        {
            er[i][j] = std::real(e1);
            ei[i][j] = std::imag(e1);
            e1 *= e0;
        }
    }

    return std::make_tuple(M, m, er, ei);
}

void fft_dif_radix2(const size_t& n, double* dr, double* di)
{
    size_t M;
    size_t *m;
    double **er;
    double **ei;

    std::tie(M, m, er, ei) = complex_factors(n, 2);

    size_t i;
    size_t j;
    size_t k;
    size_t j1;
    size_t j2;
    double tmp_r;
    double tmp_i;

    for (i = 0; i < M; ++i)
    {
        for (j = 0; j < n / (2 * m[i]); ++j)
        {
#pragma omp simd
            for (k = 0; k < m[i]; ++k)
            {
                j1 = 2 * m[i] * j + k;
                j2 = j1 + m[i];

                tmp_r = dr[j2];
                tmp_i = di[j2];

                dr[j2] = dr[j1] - tmp_r;
                di[j2] = di[j1] - tmp_i;

                dr[j1] += tmp_r;
                di[j1] += tmp_i;

                tmp_r = dr[j2];
                dr[j2] = dr[j2] * er[i][k] - di[j2] * ei[i][k];
                di[j2] = tmp_r * ei[i][k] + di[j2] * er[i][k];
            }
        }
    }
}

void fft_dif_radix4(const size_t& n, double* dr, double* di)
{
    size_t M;
    size_t *m;
    double **er;
    double **ei;

    std::tie(M, m, er, ei) = complex_factors(n, 4);

    size_t i;
    size_t j;
    size_t k;
    size_t j1;
    size_t j2;
    size_t j3;
    size_t j4;
    double tmp_r;
    double tmp_i;
    double e_r;
    double e_i;
    double a_r;
    double a_i;
    double b_r;
    double b_i;
    double c_r;
    double c_i;
    double d_r;
    double d_i;

    for (i = 0; i < M; ++i)
    {
        for (j = 0; j < n / (4 * m[i]); ++j)
        {
#pragma omp simd
            for (k = 0; k < m[i]; ++k)
            {
                j1 = 4 * m[i] * j + k;
                j2 = j1 + m[i];
                j3 = j2 + m[i];
                j4 = j3 + m[i];

                a_r = dr[j1];
                a_i = di[j1];
                b_r = dr[j2];
                b_i = di[j2];
                c_r = dr[j3];
                c_i = di[j3];
                d_r = dr[j4];
                d_i = di[j4];

                dr[j1] = a_r + b_r + c_r + d_r;
                di[j1] = a_i + b_i + c_i + d_i;

                e_r = er[i][k];
                e_i = ei[i][k];

                tmp_r = a_r + b_i - c_r - d_i;
                tmp_i = a_i - b_r - c_i + d_r;
                dr[j3] = tmp_r * e_r - tmp_i * e_i;
                di[j3] = tmp_i * e_r + tmp_r * e_i;

                tmp_i = 2 * e_r * e_i;
                e_r = e_r * e_r - e_i * e_i;
                e_i = tmp_i;

                tmp_r = a_r - b_r + c_r - d_r;
                tmp_i = a_i - b_i + c_i - d_i;
                dr[j2] = tmp_r * e_r - tmp_i * e_i;
                di[j2] = tmp_i * e_r + tmp_r * e_i;

                tmp_i = e_i * er[i][k] + e_r * ei[i][k];
                e_r = e_r * er[i][k] - e_i * ei[i][k];
                e_i = tmp_i;

                tmp_r = a_r - b_i - c_r + d_i;
                tmp_i = a_i + b_r - c_i - d_r;
                dr[j4] = tmp_r * e_r - tmp_i * e_i;
                di[j4] = tmp_i * e_r + tmp_r * e_i;
            }
        }
    }
}
