#include "workers.h"
#include "mex.h"
#include "getmba.h"
#include "dupmba.h"

/** Matlab MEX Gateway Routine for "udut_fact.f"
 *
 * function [U D] = udut_data(P)
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

    doublereal     *P;
    doublereal     *U;
    doublereal     *D;

    if (nrhs != 1)
        mexErrMsgTxt ("inputs are (P)");
    if (nlhs != 2)
        mexErrMsgTxt ("outputs are (U,D)");

    n = mxGetM (prhs[0]);

    plhs[0] = mxCreateDoubleMatrix (n, n, mxREAL);
    plhs[1] = mxCreateDoubleMatrix (n, 1, mxREAL);

    DUPMBA (P, prhs[0], n, n);
    GETMBA (U, plhs[0], n, n);
    GETMBA (D, plhs[1], n, 1);

#ifndef USE_FORTRAN
    udut_fact (n, P, U, D);
#else                           /*USE_FORTRAN */
    udut_fact__ (&n, P, U, D);
#endif                          /*USE_FORTRAN */
}
