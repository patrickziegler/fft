#include "fft.hpp"
#include <cstddef>
#include <mex.h>

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    size_t n;
    double *in_r, *in_i, *out_r, *out_i;

    if (nrhs != 1) {
        mexErrMsgTxt("Wrong number of input arguments!");
    }

    if (nlhs != 1) {
        mexErrMsgTxt("Wrong number of output arguments!");
    }

    if (mxGetM(prhs[0]) != 1) {
        mexErrMsgTxt("Input was not a row vector!");
    }

    n = static_cast<size_t>(mxGetN(prhs[0]));
    plhs[0] = mxCreateDoubleMatrix(n, 1, mxCOMPLEX);

    in_r = mxGetPr(prhs[0]);
    in_i = mxGetPi(prhs[0]);
    out_r = mxGetPr(plhs[0]);
    out_i = mxGetPi(plhs[0]);

    fft_dif_radix2(n, in_r, in_i);
    copy_bitrev(n, in_r, in_i, out_r, out_i);
}
