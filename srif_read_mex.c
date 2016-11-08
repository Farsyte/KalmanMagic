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

    doublereal     *x;
    doublereal     *P;
    doublereal     *S;
    doublereal     *Im;
    doublereal     *Iv;

    if (nrhs != 2)
        mexErrMsgTxt ("input is (Im,Iv)");
    if (nlhs != 3)
        mexErrMsgTxt ("output is (x,P,S)");

    n = mxGetNumberOfElements (prhs[1]);

    plhs[0] = mxCreateDoubleMatrix (n, 1, mxREAL);
    plhs[1] = mxCreateDoubleMatrix (n, n, mxREAL);
    plhs[2] = mxCreateDoubleMatrix (n, n, mxREAL);

    GETMBA (x, plhs[0], n, 1);
    GETMBA (P, plhs[1], n, n);
    GETMBA (S, plhs[2], n, n);

    /*
     * XXX: we can avoid the DUP,
     * * if we make sure the service code
     * * does not modify the data.
     */

    DUPMBA (Im, prhs[0], n, n);
    DUPMBA (Iv, prhs[1], n, 1);

#ifndef USE_FORTRAN
    srif_read (n, x, P, S, Im, Iv);
#else                           /*USE_FORTRAN */
    srif_read__ (&n, x, P, S, Im, Iv);
#endif                          /*USE_FORTRAN */
}
