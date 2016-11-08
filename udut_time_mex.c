#include "workers.h"
#include "mex.h"
#include "getmba.h"
#include "dupmba.h"

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
    integer         n;
    integer         m;

    doublereal     *x;
    doublereal     *U;
    doublereal     *d;
    doublereal     *F;
    doublereal     *G;

    if (nrhs != 5)
        mexErrMsgTxt ("input is (x,U,d,F,G)");
    if (nlhs != 3)
        mexErrMsgTxt ("output is (x,U,d)");

    n = mxGetM (prhs[4]);
    m = mxGetN (prhs[4]);

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);
    plhs[2] = mxDuplicateArray (prhs[2]);

    GETMBA (x, plhs[0], n, 1);
    GETMBA (U, plhs[1], n, n);
    GETMBA (d, plhs[2], n, 1);

    /*
     * XXX: make sure the worker does not destroy these.
     */
    DUPMBA (F, prhs[3], n, n);
    DUPMBA (G, prhs[4], n, m);

#ifndef USE_FORTRAN
    udut_time (n, m, x, U, d, F, G);
#else                           /*USE_FORTRAN */
    udut_time__ (&n, &m, x, U, d, F, G);
#endif                          /*USE_FORTRAN */
}
