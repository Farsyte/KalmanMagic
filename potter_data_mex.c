#include "workers.h"
#include "mex.h"
#include "getmba.h"
#include "dupmba.h"

/** Matlab MEX Gateway Routine for "potter_data.f"
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
    doublereal     *S;
    doublereal     *z;
    doublereal     *H;

    if (nrhs != 4)
        mexErrMsgTxt ("inputs are (x,S,H,z)");
    if (nlhs != 2)
        mexErrMsgTxt ("outputs are (x,S)");

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);

    n = mxGetNumberOfElements (prhs[0]);
    m = mxGetNumberOfElements (prhs[3]);

    GETMBA (x, plhs[0], n, 1);
    GETMBA (S, plhs[1], n, n);

    DUPMBA (H, prhs[2], m, n);
    DUPMBA (z, prhs[3], m, 1);

#ifndef USE_FORTRAN
    potter_data (n, m, x, S, H, z);
#else                           /*USE_FORTRAN */
    potter_data__ (&n, &m, x, S, H, z);
#endif                          /*USE_FORTRAN */
}
