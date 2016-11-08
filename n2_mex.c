#include "cfort.h"
#include "mex.h"
#include <math.h>

/** Matlab MEX Gateway Routine
 */
void
mexFunction (
    int nlhs,
    mxArray * plhs[],
    int nrhs,
    mxArray const *prhs[]
)
{
    long            n;
    long            i;
    doublereal     *p;
    double          n2;

    if (nrhs != 1)
        mexErrMsgTxt ("input is (V)");
    if (nlhs >= 2)
        mexErrMsgTxt ("output is a single scalar");

    if (!mxIsDouble (prhs[0]))
        mexErrMsgTxt ("input must be numeric.");

    if (mxIsSparse (prhs[0]))
        n = mxGetJc (prhs[0])[mxGetN (prhs[0])];
    else
        n = mxGetNumberOfElements (prhs[0]);

    n2 = 0;

    p = mxGetPr (prhs[0]);
    if (p)
        for (i = 0; i < n; ++i)
            n2 += p[i] * p[i];

    p = mxGetPi (prhs[0]);
    if (p)
        for (i = 0; i < n; ++i)
            n2 += p[i] * p[i];

    plhs[0] = mxCreateDoubleScalar (n2);

    if (plhs[0] == 0)
        mexErrMsgTxt ("Unable to construct return value");
}
