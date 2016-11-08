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

    doublereal     *Im;
    doublereal     *Iv;
    doublereal     *H;
    doublereal     *Z;

    if (nrhs != 4)
        mexErrMsgTxt ("input is (Im,Iv,H,Z)");
    if (nlhs != 2)
        mexErrMsgTxt ("output is (Im,Iv)");

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);

    n = mxGetNumberOfElements (prhs[1]);
    m = mxGetNumberOfElements (prhs[3]);

    GETMBA (Im, plhs[0], n, n);
    GETMBA (Iv, plhs[1], n, 1);

    /*
     * XXX: we can avoid the DUP,
     * * if we make sure the service code
     * * does not modify the data.
     */

    DUPMBA (H, prhs[2], m, n);
    DUPMBA (Z, prhs[3], m, 1);

#ifndef USE_FORTRAN
    srif_data (n, m, Im, Iv, H, Z);
#else                           /*USE_FORTRAN */
    srif_data__ (&n, &m, Im, Iv, H, Z);
#endif                          /*USE_FORTRAN */
}
