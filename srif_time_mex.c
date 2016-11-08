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
    doublereal     *Fi;
    doublereal     *G;
    doublereal     *Wm;
    doublereal     *Wv;

/*         [Im Iv] = srif_time_mex(Im,Iv,Fi,G,Wm,Wv);
 */

    if (nrhs != 6)
        mexErrMsgTxt ("input is (Im,Iv,Fi,G,Wm,Wv)");
    if (nlhs != 2)
        mexErrMsgTxt ("output is (Im,Iv)");

    plhs[0] = mxDuplicateArray (prhs[0]);
    plhs[1] = mxDuplicateArray (prhs[1]);

    n = mxGetM (prhs[4]);
    m = mxGetN (prhs[4]);

    GETMBA (Im, plhs[0], n, n);
    GETMBA (Iv, plhs[1], n, 1);

    /*
     * XXX: we can avoid the DUP,
     * * if we make sure the service code
     * * does not modify the data.
     */

    DUPMBA (Fi, prhs[2], n, n);
    DUPMBA (G, prhs[3], n, m);
    DUPMBA (Wm, prhs[4], n, n);
    DUPMBA (Wv, prhs[5], m, 1);

#ifndef USE_FORTRAN
    srif_time (n, m, Im, Iv, Fi, G, Wm, Wv);
#else                           /*USE_FORTRAN */
    srif_time__ (&n, &m, Im, Iv, Fi, G, Wm, Wv);
#endif                          /*USE_FORTRAN */
}
