#include "fft.hpp"
#include <complex>

struct param_radix2 {
    size_t i;
    size_t j;
    size_t k;
    size_t j1;
    size_t j2;
    double tmp_r;
    double tmp_i;
};

struct param_radix4 {
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
};

struct complex_factors {
    size_t M;
    size_t *m;
    double **er;
    double **ei;
};

void calc_complex_factors(complex_factors& c, const size_t n, const size_t radix)
{
    size_t m1, m2;
    std::complex<double> e0, e1;

    c.M = static_cast<size_t>(log2(n) / log2(radix));
    c.m = new size_t[c.M];
    c.er = new double*[c.M];
    c.ei = new double*[c.M];
    m2 = n;

    for (size_t i = 0; i < c.M; ++i)
    {
        m1 = m2;
        m2 /= radix;
        c.m[i] = m2;

        c.er[i] = new double[m2];
        c.ei[i] = new double[m2];

        e0 = std::exp(std::complex<double>(0, -2 * M_PI / m1));
        e1 = 1;

        for (size_t j = 0; j < m2; ++j)
        {
            c.er[i][j] = std::real(e1);
            c.ei[i][j] = std::imag(e1);
            e1 *= e0;
        }
    }
}

void delete_complex_factors(complex_factors& c)
{
    for (size_t i = 0; i < c.M; ++i) {
        delete[] c.er[i];
        delete[] c.ei[i];
    }

    delete[] c.er;
    delete[] c.ei;
    delete[] c.m;
}

inline void radix2(const complex_factors& c, param_radix2& p, double* dr, double* di)
{
    p.j1 = 2 * c.m[p.i] * p.j + p.k;
    p.j2 = p.j1 + c.m[p.i];

    p.tmp_r = dr[p.j2];
    p.tmp_i = di[p.j2];

    dr[p.j2] = dr[p.j1] - p.tmp_r;
    di[p.j2] = di[p.j1] - p.tmp_i;

    dr[p.j1] += p.tmp_r;
    di[p.j1] += p.tmp_i;

    p.tmp_r = dr[p.j2];
    dr[p.j2] = dr[p.j2] * c.er[p.i][p.k] - di[p.j2] * c.ei[p.i][p.k];
    di[p.j2] = p.tmp_r * c.ei[p.i][p.k] + di[p.j2] * c.er[p.i][p.k];
}

inline void radix4(const complex_factors& c, param_radix4& p, double* dr, double* di)
{
    p.j1 = 4 * c.m[p.i] * p.j + p.k;
    p.j2 = p.j1 + c.m[p.i];
    p.j3 = p.j2 + c.m[p.i];
    p.j4 = p.j3 + c.m[p.i];

    p.a_r = dr[p.j1];
    p.a_i = di[p.j1];
    p.b_r = dr[p.j2];
    p.b_i = di[p.j2];
    p.c_r = dr[p.j3];
    p.c_i = di[p.j3];
    p.d_r = dr[p.j4];
    p.d_i = di[p.j4];

    dr[p.j1] = p.a_r + p.b_r + p.c_r + p.d_r;
    di[p.j1] = p.a_i + p.b_i + p.c_i + p.d_i;

    p.e_r = c.er[p.i][p.k];
    p.e_i = c.ei[p.i][p.k];

    p.tmp_r = p.a_r + p.b_i - p.c_r - p.d_i;
    p.tmp_i = p.a_i - p.b_r - p.c_i + p.d_r;
    dr[p.j3] = p.tmp_r * p.e_r - p.tmp_i * p.e_i;
    di[p.j3] = p.tmp_i * p.e_r + p.tmp_r * p.e_i;

    p.tmp_i = 2 * p.e_r * p.e_i;
    p.e_r = p.e_r * p.e_r - p.e_i * p.e_i;
    p.e_i = p.tmp_i;

    p.tmp_r = p.a_r - p.b_r + p.c_r - p.d_r;
    p.tmp_i = p.a_i - p.b_i + p.c_i - p.d_i;
    dr[p.j2] = p.tmp_r * p.e_r - p.tmp_i * p.e_i;
    di[p.j2] = p.tmp_i * p.e_r + p.tmp_r * p.e_i;

    p.tmp_i = p.e_i * c.er[p.i][p.k] + p.e_r * c.ei[p.i][p.k];
    p.e_r = p.e_r * c.er[p.i][p.k] - p.e_i * c.ei[p.i][p.k];
    p.e_i = p.tmp_i;

    p.tmp_r = p.a_r - p.b_i - p.c_r + p.d_i;
    p.tmp_i = p.a_i + p.b_r - p.c_i - p.d_r;
    dr[p.j4] = p.tmp_r * p.e_r - p.tmp_i * p.e_i;
    di[p.j4] = p.tmp_i * p.e_r + p.tmp_r * p.e_i;
}

void fft_dif_radix2(const size_t n, double* dr, double* di)
{
    param_radix2 p;
    complex_factors c;

    p.i = 0;
    p.j = 0;
    p.k = 0;

    calc_complex_factors(c, n, 2);

    for (p.i = 0; p.i < c.M; ++p.i)
    {
        for (p.j = 0; p.j < n / (2 * c.m[p.i]); ++p.j)
        {
            for (p.k = 0; p.k < c.m[p.i]; ++p.k)
            {
                radix2(c, p, dr, di);
            }
        }
    }

    delete_complex_factors(c);
}

void fft_dif_radix4(const size_t n, double* dr, double* di)
{
    param_radix4 p;
    complex_factors c;

    p.i = 0;
    p.j = 0;
    p.k = 0;

    calc_complex_factors(c, n, 4);

    for (p.i = 0; p.i < c.M; ++p.i)
    {
        for (p.j = 0; p.j < n / (4 * c.m[p.i]); ++p.j)
        {
            for (p.k = 0; p.k < c.m[p.i]; ++p.k)
            {
                radix4(c, p, dr, di);
            }
        }
    }

    delete_complex_factors(c);
}

void copy_bitrev(size_t n, double* in_r, double* in_i, double* out_r, double* out_i)
{
    size_t tmp, m, n2 = n / 2;
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
