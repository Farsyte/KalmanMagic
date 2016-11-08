#include "workers.h"
#include "mex.h"
#include "getmba.h"

/** Matlab MEX Gateway Routine for "hh_tri.f"
 */
void
mexFunction (
    int nlhs,
    mxArray * plhs[],
    int nrhs,
    mxArray const *prhs[]
)
{
    integer         n;
    integer         m;
    doublereal     *A;
    doublereal     *D;
    integer         d;

    if (nrhs != 2)
        mexErrMsgTxt ("inputs are (A,d)");
    if (nlhs != 1)
        mexErrMsgTxt ("outputs are (A)");

    plhs[0] = mxDuplicateArray (prhs[0]);

    n = mxGetN (plhs[0]);
    m = mxGetM (plhs[0]);

    GETMBA (A, plhs[0], m, n);
    GETMBA (D, prhs[1], 1, 1);

    d = (int)(*D + 0.5);

#ifndef USE_FORTRAN
    hh_tri (m, n, A, d);
#else                           /*USE_FORTRAN */
    hh_tri__ (&m, &n, A, &d);
#endif                          /*USE_FORTRAN */
}
