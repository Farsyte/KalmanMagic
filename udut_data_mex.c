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
    doublereal     *D;
    doublereal     *z;
    doublereal     *H;
    doublereal     *R;

    setbuf (stderr, 0);

    if (nrhs != 6)
        mexErrMsgTxt ("input is (X,U,D,H,R,Z)");
    if (nlhs != 3)
        mexErrMsgTxt ("output is (X,U,D)");

    n = mxGetNumberOfElements (prhs[0]);
    m = mxGetNumberOfElements (prhs[5]);

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);
    plhs[2] = mxDuplicateArray (prhs[2]);

    GETMBA (x, plhs[0], n, 1);
    GETMBA (U, plhs[1], n, n);
    GETMBA (D, plhs[2], n, 1);

    /*
     * XXX: make sure the worker does not destroy these.
     */
    DUPMBA (H, prhs[3], m, n);
    DUPMBA (R, prhs[4], m, 1);
    DUPMBA (z, prhs[5], m, 1);

#ifndef USE_FORTRAN
    udut_data (n, m, x, U, D, H, R, z);
#else                           /*USE_FORTRAN */
    udut_data__ (&n, &m, x, U, D, H, R, z);
#endif                          /*USE_FORTRAN */
}
